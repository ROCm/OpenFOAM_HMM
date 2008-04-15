/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

// Foam header files
#include "IFstream.H"
#include "OSspecific.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "NameServer.H"
#include "LogManager.H"
#include "LogEntry.H"
#include "Paths.H"
#include "Pstream.H"

// Declare all the namespaces used by FoamX
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::NameServer::NameServer()
:
    connected_(false)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::NameServer::NameServer(CORBA::ORB_ptr pOrb)
:
    connected_(false)
{
    connect(pOrb);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::NameServer::~NameServer()
{
    disconnect();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::NameServer::connect(CORBA::ORB_ptr pOrb)
{
    static const char* functionName =
        "FoamX::NameServer::connect(CORBA::ORB_ptr pOrb)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (connected_ || !CORBA::is_nil(rootContext_))
        {
            throw FoamXError
            (
                E_FAIL,
                "Name server already connected.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Try to connect to the name server using the
        // resolve_initial_references call provided by the ORB.
        // Different ORB vendors implement the way in which their ORB obtains
        // initial references in different ways (way to go OMG boys, an
        // inspired design decision!).
        // Therefore, we cannot rely on this call working and may have to fall
        // back on the tried and trusted method of retrieving the name service
        // IOR from a file written to a common location. Since the MICO name 
        // server reference cannot be obtained from JAVA in any other way,
        // we have to write this file anyway. MICO C++ servers can specify the
        // location of the name server by using the
        // "-ORBNamingAddr inet:hostname:portno" command line argument
        // in which case the resolve_initial_references will work.
        CORBA::Object_var initServ;
        try
        {
            // Obtain a reference to the root context of the Name service
            // through the given Orb.
            // Attach to reference.
            initServ = pOrb->resolve_initial_references("NameService");

            log << "Obtained NameService reference from Orb command line"
                << " arguments" << endl;
        }
        catch (...)
        {
            // If the resolve_initial_references call failed, try and open 
            // the name service reference written to the FoamX system directory.
            if (CORBA::is_nil(initServ))
            {
                // See if the name service reference file exists.
                fileName nsRefFile = Paths::tmp/"ns.ref";
                if (!exists(nsRefFile))
                {
                    throw FoamXError
                    (
                        E_FAIL,
                        "Failed to resolve NameService root context reference.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                // Extract IOR from the file.
                log << "Reading name service reference file '"
                    << nsRefFile << "'." << endl;

                word ior((IFstream(nsRefFile)()));
                initServ = pOrb->string_to_object(ior.c_str());
            }
        }

        // If we still haven't got a name server reference, throw an exception.
        if (CORBA::is_nil(initServ))
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to resolve NameService root context reference.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // call to a CosNaming::NamingContext object.
        // Attach to reference.
        rootContext_ = CosNaming::NamingContext::_narrow(initServ);
        if (CORBA::is_nil(rootContext_))
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to narrow NameService root context reference.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Set connected flag.
        connected_ = true;
    }
    catch (CORBA::ORB::InvalidName& ex)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught CORBA::ORB::InvalidName exception in connect : "
            "Service required is invalid [does not exist].",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CORBA::COMM_FAILURE& ex)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught CORBA::COMM_FAILURE exception in connect : "
            "Unable to contact name server.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CORBA::NO_RESOURCES& ex)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught CORBA::NO_RESOURCES exception in connect : "
            "Unable to contact name server.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CORBA::SystemException& ex)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught CORBA::SystemException exception in connect : "
            "Unable to contact name server.",
            functionName, 
            __FILE__, __LINE__
        );
    }
    catch (...)
    {
        throw FoamXError
        (
            E_UNEXPECTED,
            "Unexpected error",
            functionName,
            __FILE__, __LINE__
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::NameServer::disconnect()
{
    static const char* functionName =
        "FoamX::NameServer::disconnect()";

    LogEntry log(functionName, __FILE__, __LINE__);

    if (connected_)
    {
        // Dicsonnect from the root context and reset the connected flag.
        rootContext_ = CosNaming::NamingContext::_nil();
        connected_   = false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::NameServer::createContext
(
    const fileName& path,
    bool failAlreadyBound
)
{
    static const char* functionName =
        "FoamX::NameServer::createContext"
        "(const fileName& path, bool failAlreadyBound)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (!connected_ || CORBA::is_nil(rootContext_))
        {
            throw FoamXError
            (
                E_FAIL,
                "Name server not connected.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Construct a CosNaming::Name from the given path string.
        CosNaming::Name cosName;
        createNameFromString(path, cosName);

        // Call bind_new_context on the root context.
        // Attach to reference.
        CosNaming::NamingContext_var context =
            rootContext_->bind_new_context(cosName);
    }
    catch (CosNaming::NamingContext::AlreadyBound)
    {
        // Raise an exception only if required.
        if (failAlreadyBound)
        {
            throw FoamXError
            (
                E_FAIL,
                "Caught a CosNaming::NamingContext::AlreadyBound exception.",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    catch (CosNaming::NamingContext::NotFound& ex)
    {
        string msg = "CosNaming::NamingContext::NotFound Exception : ";

        if (ex.why == CosNaming::NamingContext::missing_node)
        {
            msg += "Missing Node.";
        }
        else if (ex.why == CosNaming::NamingContext::not_context)
        {
            msg += "Not Context.";
        }
        else if (ex.why == CosNaming::NamingContext::not_object)
        {
            msg += "Not Object.";
        }
        throw FoamXError
        (
            E_FAIL,
            msg,
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::InvalidName)
    {
         throw FoamXError
         (
             E_FAIL,
             "Caught a CosNaming::NamingContext::InvalidName exception.",
             functionName,
             __FILE__, __LINE__
         );
    }
    catch (CosNaming::NamingContext::NotEmpty)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::NotEmpty exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::CannotProceed)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::CannotProceed exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CORBA::COMM_FAILURE)
    {
         throw FoamXError
         (
             E_FAIL,
             "Caught a CORBA::COMM_FAILURE exception.",
             functionName,
             __FILE__, __LINE__
         );
    }
    catch (...)
    {
         throw FoamXError
         (
             E_FAIL,
             "Unexpected error.",
             functionName,
             __FILE__, __LINE__
         );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::NameServer::createContexts
(
    const fileName& path,
    bool failAlreadyBound
)
{
    static const char* functionName =
        "FoamX::NameServer::createContext"
        "(const fileName& path, bool failAlreadyBound)";

    LogEntry log(functionName, __FILE__, __LINE__);


    // Pass errors through

    wordList components = path.components();

    fileName dir = "";

    forAll(components, compI)
    {
        if (compI == 0)
        {
            dir = components[compI];
        }
        else
        {
            dir = dir/components[compI];
        }
        createContext(dir, failAlreadyBound);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::NameServer::removeContext(const fileName& path)
{
    static const char* functionName =
        "FoamX::NameServer::removeContext(const fileName& path)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (!connected_ || CORBA::is_nil(rootContext_))
        {
            throw FoamXError
            (
                E_FAIL,
                "Name server not connected.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Construct a CosNaming::Name from the given path string.
        CosNaming::Name cosName;
        createNameFromString(path, cosName);

        // Get the specified naming context object.
        CORBA::Object_var obj = rootContext_->resolve(cosName);
        CosNaming::NamingContext_var context =
            CosNaming::NamingContext::_narrow(obj);

        if (CORBA::is_nil(context))
        {
            throw FoamXError
            (
                E_FAIL,
                "Not a naming context.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Destroy and unbind the specified naming context.
        context->destroy();
        rootContext_->unbind(cosName);
    }
    catch (CosNaming::NamingContext::NotFound& ex)
    {
        string msg = "CosNaming::NamingContext::NotFound Exception : ";

        if (ex.why == CosNaming::NamingContext::missing_node)
        {
            msg += "Missing node.";
        }
        else if (ex.why == CosNaming::NamingContext::not_context)
        {
            msg += "Not Context.";
        }
        else if (ex.why == CosNaming::NamingContext::not_object)
        {
            msg += "Not Object.";
        }
        throw FoamXError
        (
            E_FAIL,
            msg,
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::InvalidName)
    {
         throw FoamXError
         (
             E_FAIL,
             "Caught a CosNaming::NamingContext::InvalidName exception.",
             functionName,
             __FILE__, __LINE__
         );
    }
    catch (CosNaming::NamingContext::AlreadyBound)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::AlreadyBound exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::NotEmpty)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::NotEmpty exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::CannotProceed)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::CannotProceed exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CORBA::COMM_FAILURE)
    {
         throw FoamXError
         (
             E_FAIL,
             "Caught a CORBA::COMM_FAILURE exception.",
             functionName,
             __FILE__, __LINE__
         );
    }
    catch (...)
    {
         throw FoamXError
         (
             E_FAIL,
             "Unexpected error.",
             functionName,
             __FILE__, __LINE__
         );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::NameServer::removeContexts
(
    const fileName& root,
    const fileName& path
)
{
    static const char* functionName =
        "FoamX::NameServer::removeContexts"
        "(const fileName& root, const fileName& path)";

    LogEntry log(functionName, __FILE__, __LINE__);

    fileName dir(path);

    while((dir != "/") && (dir != "") && (dir != root))
    {
        removeContext(dir);
        dir = dir.path();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::NameServer::bindObject
(
    const fileName& name,
    CORBA::Object_ptr pObject,
    bool rebind
)
{
    static const char* functionName =
        "FoamX::NameServer::bindObject"
        "(const fileName& name, CORBA::Object_ptr pObject, bool rebind)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (!connected_ || CORBA::is_nil(rootContext_))
        {
            throw FoamXError
            (
                E_FAIL,
                "Name server not connected.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Construct a CosNaming::Name from the given path string.
        CosNaming::Name cosName;
        createNameFromString(name, cosName);

        // Bind this object to the specified name.
        try
        {
            rootContext_->bind(cosName, pObject);
        }
        catch (CosNaming::NamingContext::AlreadyBound& ex)
        {
            if (!rebind) throw ex;
            // Note: Using rebind() will overwrite any Object previously bound
            //       to the host context with obj.
            //       Alternatively, bind() can be used, which will raise a
            //       CosNaming::NamingContext::AlreadyBound exception if the
            //       name supplied is already bound to an object.
            // Amendment: When using OrbixNames, it is necessary to first try
            // bind and then rebind, as rebind on it's own will throw a
            // NotFoundexception if the Name has not already been bound.
            //[This is incorrect behaviour - it should just bind].
            rootContext_->rebind(cosName, pObject);
        }
    }
    catch (CosNaming::NamingContext::NotFound& ex)
    {
        string msg = "CosNaming::NamingContext::NotFound Exception : ";

        if (ex.why == CosNaming::NamingContext::missing_node)
        {
            msg += "Missing node.";
        }
        else if (ex.why == CosNaming::NamingContext::not_context)
        {
            msg += "Not Context.";
        }
        else if (ex.why == CosNaming::NamingContext::not_object)
        {
            msg += "Not Object.";
        }
        throw FoamXError
        (
            E_FAIL,
            msg,
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::InvalidName)
    {
         throw FoamXError
         (
             E_FAIL,
             "Caught a CosNaming::NamingContext::InvalidName exception.",
             functionName,
             __FILE__, __LINE__
         );
    }
    catch (CosNaming::NamingContext::AlreadyBound)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::AlreadyBound exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::NotEmpty)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::NotEmpty exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::CannotProceed)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::CannotProceed exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CORBA::COMM_FAILURE)
    {
         throw FoamXError
         (
             E_FAIL,
             "Caught a CORBA::COMM_FAILURE exception.",
             functionName,
             __FILE__, __LINE__
         );
    }
    catch (...)
    {
         throw FoamXError
         (
             E_FAIL,
             "Unexpected error.",
             functionName,
             __FILE__, __LINE__
         );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::NameServer::unbindObject(const fileName& name)
{
    static const char* functionName =
        "FoamX::NameServer::unbindObject(const fileName& name)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (!connected_ || CORBA::is_nil(rootContext_))
        {
            throw FoamXError
            (
                E_FAIL,
                "Name server not connected.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Construct a CosNaming::Name from the given path string.
        CosNaming::Name cosName;
        createNameFromString(name, cosName);

        // List parent context and check that binding type is "nobject"
        // before allowing unbind.
        CosNaming::NamingContext_var context;
        int len = cosName.length();
        CORBA::String_var id   = cosName[len-1].id;
        CORBA::String_var kind = cosName[len-1].kind;

        // Check that the parent is a naming context.
        if (len> 1)
        {
            cosName.length(len - 1);
            CORBA::Object_var obj = rootContext_->resolve(cosName);
            context = CosNaming::NamingContext::_narrow(obj);
            if (CORBA::is_nil(context))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Parent is not a naming context.",
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }
        else
        {
            context = rootContext_;
        }

        // Loop over all bindings in the parent context.
        CosNaming::BindingIterator_var bi;
        CosNaming::BindingList_var     bl;
        CosNaming::Binding_var         b;

        bool bFound = false;
        context->list(0, bl, bi);
        while (bi->next_one(b))
        {
            if
            (
                strcmp(b->binding_name[0].id, id   ) == 0
             && strcmp(b->binding_name[0].kind, kind) == 0
            )
            {
                if (b->binding_type == CosNaming::ncontext)
                {
                    throw FoamXError
                    (
                        E_FAIL,
                        "Can't unbind a naming context.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                // Unbind the object.
                context->unbind(b->binding_name);
                bFound = true;
                break;
            }
        }
        bi->destroy();          //?

        if (!bFound)
        {
            throw FoamXError
            (
                E_FAIL,
                "Object name could not be found.",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    catch (CosNaming::NamingContext::NotFound& ex)
    {
        string msg = "CosNaming::NamingContext::NotFound Exception : ";

        if (ex.why == CosNaming::NamingContext::missing_node)
        {
            msg += "Missing node.";
        }
        else if (ex.why == CosNaming::NamingContext::not_context)
        {
            msg += "Not Context.";
        }
        else if (ex.why == CosNaming::NamingContext::not_object)
        {
            msg += "Not Object.";
        }
        throw FoamXError
        (
            E_FAIL,
            msg,
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::InvalidName)
    {
         throw FoamXError
         (
             E_FAIL,
             "Caught a CosNaming::NamingContext::InvalidName exception.",
             functionName,
             __FILE__, __LINE__
         );
    }
    catch (CosNaming::NamingContext::AlreadyBound)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::AlreadyBound exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::NotEmpty)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::NotEmpty exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::CannotProceed)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::CannotProceed exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CORBA::COMM_FAILURE)
    {
         throw FoamXError
         (
             E_FAIL,
             "Caught a CORBA::COMM_FAILURE exception.",
             functionName,
             __FILE__, __LINE__
         );
    }
    catch (...)
    {
         throw FoamXError
         (
             E_FAIL,
             "Unexpected error.",
             functionName,
             __FILE__, __LINE__
         );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool FoamX::NameServer::isObjectBound(const fileName& name)
{
    static const char* functionName =
        "FoamX::NameServer::isObjectBound(const fileName& name)";

    LogEntry log(functionName, __FILE__, __LINE__);

    bool objectBound = false;

    try
    {
        if (!connected_ || CORBA::is_nil(rootContext_))
        {
            throw FoamXError
            (
                E_FAIL,
                "Name server not connected.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Construct a CosNaming::Name from the given path string.
        CosNaming::Name cosName;
        createNameFromString(name, cosName);

        // List parent context and check that binding type is "nobject".
        CosNaming::NamingContext_var context;
        int len = cosName.length();
        CORBA::String_var id   = cosName[len-1].id;
        CORBA::String_var kind = cosName[len-1].kind;

        // Check that the parent is a naming context.
        if (len> 1)
        {
            cosName.length(len - 1);

            // Attach to reference.
            CORBA::Object_var obj = rootContext_->resolve(cosName);

            // Attach to reference.
            context = CosNaming::NamingContext::_narrow(obj);
            if (CORBA::is_nil(context))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Parent is not a naming context.",
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }
        else
        {
            context = rootContext_;
        }

        // Loop over all bindings in the parent context.
        CosNaming::BindingIterator_var bi;
        CosNaming::BindingList_var     bl;
        CosNaming::Binding_var         b;

        context->list(0, bl, bi);

        if (bi)
        {
            while (bi->next_one(b))
            {
                if
                (
                    strcmp(b->binding_name[0].id, id   ) == 0
                 && strcmp(b->binding_name[0].kind, kind) == 0
                )
                {
                    objectBound = b->binding_type == CosNaming::nobject;
                    break;
                }
            }
            bi->destroy();
        }
    }
    catch (CosNaming::NamingContext::NotFound& ex)
    {
        string msg = "CosNaming::NamingContext::NotFound Exception : ";

        if (ex.why == CosNaming::NamingContext::missing_node)
        {
            // msg += "Missing node.";
            // Missing node is no longer treated as an exception,
            // it simply implies the object is not bound
            return false;
        } 
        else if (ex.why == CosNaming::NamingContext::not_context)
        {
            msg += "Not Context.";
        } 
        else if (ex.why == CosNaming::NamingContext::not_object)
        {
            msg += "Not Object.";
        }
        throw FoamXError
        (
            E_FAIL,
            msg,
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::InvalidName)
    {
         throw FoamXError
         (
             E_FAIL,
             "Caught a CosNaming::NamingContext::InvalidName exception.",
             functionName,
             __FILE__, __LINE__
         );
    }
    catch (CosNaming::NamingContext::AlreadyBound)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::AlreadyBound exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::NotEmpty)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::NotEmpty exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::CannotProceed)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::CannotProceed exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CORBA::COMM_FAILURE)
    {
         throw FoamXError
         (
             E_FAIL,
             "Caught a CORBA::COMM_FAILURE exception.",
             functionName,
             __FILE__, __LINE__
         );
    }
    catch (...)
    {
         /* This is commented out because if you open a FoamXCaseBrowser
            close it and then reopen it this error is thrown.  Really we need
            to find out what exception is actually being thrown and catch it
         throw FoamXError
         (
             E_FAIL,
             "Unexpected error.",
             functionName,
             __FILE__, __LINE__
         );
         */

         return false;
    }

    return objectBound;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Object_ptr FoamX::NameServer::resolve(const fileName& name)
{
    static const char* functionName =
        "FoamX::NameServer::resolve(const fileName& name)";

    LogEntry log(functionName, __FILE__, __LINE__);

    CORBA::Object_ptr pObject = CORBA::Object::_nil();

    try
    {
        if (!connected_ || CORBA::is_nil(rootContext_))
        {
            throw FoamXError
            (
                E_FAIL,
                "Name server not connected.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Construct a CosNaming::Name from the given path string.
        CosNaming::Name cosName;
        createNameFromString(name, cosName);

        // Resolve the name and get the object reference.
        // Attach to reference.
        CORBA::Object_var obj = rootContext_->resolve(cosName);

        // Duplicate the reference and return.
        pObject = CORBA::Object::_duplicate(obj);
    }
    catch (CosNaming::NamingContext::NotFound& ex)
    {
        string msg = "CosNaming::NamingContext::NotFound Exception : ";

        if (ex.why == CosNaming::NamingContext::missing_node)
        {
            msg += "Missing node.";
        }
        else if (ex.why == CosNaming::NamingContext::not_context)
        {
            msg += "Not Context.";
        }
        else if (ex.why == CosNaming::NamingContext::not_object)
        {
            msg += "Not Object.";
        }
        throw FoamXError
        (
            E_FAIL,
            msg,
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::InvalidName)
    {
         throw FoamXError
         (
             E_FAIL,
             "Caught a CosNaming::NamingContext::InvalidName exception.",
             functionName,
             __FILE__, __LINE__
         );
    }
    catch (CosNaming::NamingContext::AlreadyBound)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::AlreadyBound exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::NotEmpty)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::NotEmpty exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CosNaming::NamingContext::CannotProceed)
    {
        throw FoamXError
        (
            E_FAIL,
            "Caught a CosNaming::NamingContext::CannotProceed exception.",
            functionName,
            __FILE__, __LINE__
        );
    }
    catch (CORBA::COMM_FAILURE)
    {
         throw FoamXError
         (
             E_FAIL,
             "Caught a CORBA::COMM_FAILURE exception.",
             functionName,
             __FILE__, __LINE__
         );
    }
    catch (...)
    {
         throw FoamXError
         (
             E_FAIL,
             "Unexpected error.",
             functionName,
             __FILE__, __LINE__
         );
    }

    return pObject;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// "Borrowed" from the MICO NameServerClient class.

void FoamX::NameServer::createNameFromString
(
    const std::string& str,
    CosNaming::Name& name
)
{
    size_t pos = 0, p = 0;
    CORBA::ULong num = 0;
    if (str[0] == '/')
    {
        pos = 1;
    }

    do
    {
        p = str.find('/', pos);
        num++;

        std::string sub;

        if (CORBA::Long(p) <0)
        {
            sub = str.substr(pos);
        }
        else
        {
            sub = str.substr(pos, p - pos);
        }

        pos = p + 1;
        name.length(num);
        name[num - 1].id   = CORBA::string_dup(sub.c_str());
        name[num - 1].kind = CORBA::string_dup("");

    } while (CORBA::Long(p) >= 0);
}


// ************************************************************************* //

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

#include <unistd.h>

// Foam header files.
#include "long.H"
#include "dictionary.H"
#include "wordList.H"
#include "stringList.H"
#include "fileNameList.H"
#include "IFstream.H"
#include "OStringStream.H"
#include "instantList.H"
#include "IOPtrList.H"
#include "OSspecific.H"

// FoamX header files.
#include "FoamXErrors.H"
#include "ICaseServerImpl.H"
#include "DictionaryWriter.H"
#include "RootDictionary.H"
#include "IGeometricFieldImpl.H"
#include "LogEntry.H"
#include "NameServer.H"
#include "PatchProperties.H"
#include "Orb.H"
#include "Paths.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Define word for IOPtrList<entry>.
template<>
const Foam::word Foam::IOPtrList<Foam::entry>::typeName("polyBoundaryMesh");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Time& FoamX::ICaseServerImpl::time()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::time()";

    if (!dbPtr_)
    {
        throw FoamXError
        (
            E_FAIL,
            "database not yet constructed",
            functionName,
            __FILE__, __LINE__
        );
    }
    return *dbPtr_;
}


void FoamX::ICaseServerImpl::saveControlDict()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::saveControlDict()";

    LogEntry log(functionName, __FILE__, __LINE__);

    //----------------------------------------------------------------------
    //  Make sure the control dictionary is saved.
    IDictionaryEntry_var controlDict;      // Auto Release.
    getDictionary(Time::controlDictName.c_str(), false, controlDict.out());
    if (CORBA::is_nil(controlDict))
    {
        throw FoamXError
        (
            E_FAIL,
            "Failed to open control dictionary entry object",
            functionName,
            __FILE__, __LINE__
        );
    }

    // Make sure the control dictionary is saved.
    controlDict->save();
}


void FoamX::ICaseServerImpl::writePatchData()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::writePatchData()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Create the polyMesh directory if required.
        fileName polyMeshDir = rootDir_/caseName_/"constant"/"polyMesh";
        if (!dir(polyMeshDir) && !mkDir(polyMeshDir))
        {
            throw FoamXError
            (
                E_FAIL,
                "Mesh directory '" +  polyMeshDir + "' could not be created",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Write boundary file (fixed form).
        DictionaryWriter dict
        (
            rootDir_,
            caseName_,
            "constant/polyMesh",
            "boundary"
        );
        dict.writeHeader
        (
            "FoamX Mesh Description File",
            "polyBoundaryMesh"
        );

        // Start list.
        dict.startList(patchMap_.size());

        // Write patch boundary information.
        for
        (
            Dictionary<PatchProperties>::iterator iter = patchMap_.begin();
            iter != patchMap_.end();
            ++iter
        )
        {
            dict.writeEndl();
            iter().save(dict);
        }

        // End list.
        dict.endList();

        dict.writeEndl();
        dict.writeEndBar();
    }
    CATCH_ALL(functionName);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// patchPhysicalType is the name of the boundary type name, not the definition key.

void FoamX::ICaseServerImpl::addPatch
(
    const char* patchName,
    const char* patchPhysicalType
)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::addPatch"
        "(const char* patchName, const char* patchPhysicalType)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Patch = " << patchName
            << ", Boundary (Physical) Type = " << patchPhysicalType << endl;

        // Create a new patch descriptor object for this patch.
        PatchProperties* patchProps = new PatchProperties(patchName);

        if (patchProps == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to create Patch Properties object",
                functionName,
                __FILE__, __LINE__
            );
        }

        IPatchPhysicalTypeDescriptor_var boundaryDescriptor;
        app_->getPatchPhysicalType(patchPhysicalType, boundaryDescriptor.out());
        patchProps->patchType(boundaryDescriptor->patchType());

        patchProps->physicalType(patchPhysicalType);

        // Append to patch descriptor list.
        addPatch(patchProps);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool FoamX::ICaseServerImpl::fieldsMatchPatchPhysicalType(const word& patchName)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::fieldsMatchPatchPhysicalType"
        "(const char* patchName)";

    // Get properties read from boundary file.
    PatchProperties* patchProps = patchMap_.lookup(patchName);
    
    Info<< "Trying to match patch " << patchName << " physical type "
        << patchProps->physicalType() << endl;

    // Get the specified boundary type descriptor object (ultimately from the
    // application class.
    IPatchPhysicalTypeDescriptor_var boundaryDescriptor;
    app_->getPatchPhysicalType
    (
        patchProps->physicalType().c_str(), 
        boundaryDescriptor.out()
    );

    // Get the types for individual fields
    StringPairList_var patchFieldTypes =
        boundaryDescriptor->patchFieldTypes();

    // Try to see if all fields specified have the correct type
    for (unsigned int i = 0; i <patchFieldTypes->length(); i++)
    {
        word fieldName(patchFieldTypes[i].name);

        Info<< "    name:" << fieldName
            << "  type:" << patchFieldTypes[i].value << endl;


        if (fieldValueMap_.found(fieldName))
        {
            // Convert between the patch field type key and the Foam name.

            // Get the PatchFieldDescriptor object.
            ITypeDescriptor_var pfDescriptor;
            foamProperties_->getPatchFieldType
            (
                patchFieldTypes[i].value,
                pfDescriptor.out()
            );
            CORBA::String_var pfName = pfDescriptor->name();

            IGeometricFieldImpl* pFieldValue = fieldValueMap_[fieldName];

            IDictionaryEntry_var pfParams;
            pFieldValue->getPatchFieldParameters(patchName.c_str(), pfParams);

            DictionaryEntryList* subElemsPtr = pfParams->subElements();

            for
            (
                unsigned int elemI = 0;
                elemI < subElemsPtr->length();
                elemI++
            )
            {
                FoamXServer::FoamXAny* subValPtr =
                    (*subElemsPtr)[elemI]->value();

                void* nasty = &(subValPtr->value);

                Info<< "    " << elemI << "type:" << subValPtr->type
                    << "  value:" << long(nasty) << endl;
            }
        }
        else
        {
            WarningIn(functionName)
                << "Did not find field " << fieldName
                << " specified in physical patch type "
                << patchProps->physicalType()
                << endl;
        }
    }
    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

FoamX::ICaseServerImpl::ICaseServerImpl
(
    Orb& orb,
    const fileName& rootDir,
    const fileName& caseName,
    const word& mode,
    const word& app
)
:
    rootDir_(rootDir),
    caseName_(caseName),
    hostContext_(hostName()),
    userContext_(userName()),
    objectName_(hostContext_/userContext_/string::validate<word>(rootDir_)/caseName_),
    caseDictName_(rootDir_/caseName_/"system/controlDict"),
    appName_(app),
    orb_(orb),
    dbPtr_(NULL),
    caseBrowser_(NULL),
    foamProperties_(NULL),
    app_(NULL),
    managed_(false),
    procControl_(dotFoam("apps/FoamX/FoamX.cfg")),
    pid_(-1)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::ICaseServerImpl"
        "(Orb& orb, const fileName& rootDir, const fileName& caseName, "
        "const word& mode, const word& app)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Connect to name server.
        NameServer fxNameServer(orb_.orbPtr_);

        // Case browser should have registered itself with the name server under
        // a "hostName/"FoamXCaseBrowser"" key.
        fileName caseBrowserKey = hostName()/"FoamXCaseBrowser";

        if (fxNameServer.isObjectBound(caseBrowserKey))
        {
            // Return server reference.
            caseBrowser_ = 
                fxNameServer.resolve<FoamXServer::CaseBrowser::ICaseBrowser>
                (caseBrowserKey);
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Server '" + caseBrowserKey + "' not found on name server. "
                "getServerReference call timed out",
                functionName,
                __FILE__, __LINE__
            );
        }

        foamProperties_ = caseBrowser_->foamProperties();


        //----------------------------------------------------------------------
        // See if we are opening an existing case or creating/importing
        // a new case.
        if (mode == "open")
        {
            // Control dictionary must exist.
            if (!exists(caseDictName_))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Case control dictionary '" + caseDictName_ + "' not found",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Open the case configuration dictionary and get the
            // application class name.
            dictionary caseDict((IFstream(caseDictName_)()));

            if (caseDict.found("application"))
            {
                caseDict.lookup("application")>> appName_;
            }
            else if (caseDict.found("applicationClass"))
            {
                caseDict.lookup("applicationClass")>> appName_;
            }
            else
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Application class entry not found in '" + caseDictName_
                  + "'",
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }
        else if (mode == "import")
        {
            // Control dictionary must exist.
            if (!exists(caseDictName_))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Case control dictionary '" + caseDictName_
                  + "' not found",
                    functionName,
                    __FILE__, __LINE__
               );
            }
        
            // Open the case configuration dictionary and get the
            // application class name.
            dictionary caseDict((IFstream(caseDictName_)()));
        
            // If not found, use the command line argument.
            if (caseDict.found("application"))
            {
                caseDict.lookup("application") >> appName_;
            }
        }
        else if (mode == "create")
        {
            // Control dictionary must not exist.
            if (exists(caseDictName_))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Case control dictionary '" + caseDictName_
                  + "' exists. Cannot create new case",
                    functionName,
                    __FILE__, __LINE__
                );
            }
            // Create the case directory if required.
            fileName caseDir = rootDir_/caseName_;
            if (!dir(caseDir) && !mkDir(caseDir))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Case directory '" +  caseDir + "' could not be created",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create the constant directory if required.
            fileName constDir = caseDir/"constant";
            if (!dir(constDir ) && !mkDir(constDir ) )
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Case directory '" +  constDir + "' could not be created",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create the polyMesh directory if required.
            fileName meshDir = constDir/"polyMesh";
            if (!dir(meshDir) && !mkDir(meshDir))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Case directory '" +  meshDir + "' could not be created",
                    functionName,
                    __FILE__, __LINE__ 
                );
            }
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Invalid argument '" +  mode + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        //----------------------------------------------------------------------
        // Create the Application and FoamProperties objects.
        if (appName_.size() == 0)
        {
            throw FoamXError
            (
                E_FAIL,
                "Invalid application class name '" + appName_ + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        foamProperties_->getApplication(appName_.c_str(), app_);

        //----------------------------------------------------------------------
        // Create the FieldValues objects.
        StringList_var fieldNames = app_->fields();

        for (unsigned int i = 0; i <fieldNames->length(); i++)
        {
            word fieldName(fieldNames[i]);

            // Get reference to field descriptor.
            IGeometricFieldDescriptor_var fieldDescriptorRef; // Auto-release.
            app_->getField(fieldName.c_str(), fieldDescriptorRef.out());

            // Create default FieldValues object for this field.
            IGeometricFieldImpl* pFieldValues = new IGeometricFieldImpl
            (
                fieldDescriptorRef,
                foamProperties_
            );

            if (pFieldValues == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Couldn't create field value object",
                    functionName,
                    __FILE__, __LINE__
                );
            }
            fieldValueMap_.insert(fieldName, pFieldValues);
        }

        //----------------------------------------------------------------------
        // Register object with the naming service. Will throw an exception
        // if any errors occurs.

        // Do not fail if context already exists.
        fxNameServer.createContext(hostContext_, false);

        // Do not fail if context already exists.
        fxNameServer.createContext(hostContext_/userContext_, false);

        // Do not fail if context already exists.
        fxNameServer.createContext
        (
            hostContext_/userContext_/string::validate<word>(rootDir_),
            false
        );

        // Bind object under key hostName/userName/caseroot/casename.
        FoamXServer::CaseServer::ICaseServer_var caseRef = _this();

        // Throw exception if already bound.
        fxNameServer.bindObject(objectName_, caseRef, false);

        // Decrement the reference count of the object implementation, so
        // that it will be properly cleaned up when the POA has determined
        // that it is no longer needed.
        _remove_ref();

        if (mode == "import" || mode == "create")
        {
            // Make valid control dict so we can instantiate a database
            saveControlDict();

            dbPtr_ = new Time(Time::controlDictName, rootDir_, caseName_);

            // Make valid set of files
            save();

            caseBrowser_->addCase
            (
                rootDir_.c_str(),
                rootDir_.c_str(),
                caseName_.c_str(),
                appName_.c_str()
            );
        }
        else
        {
            dbPtr_ = new Time(Time::controlDictName, rootDir_, caseName_);

            caseBrowser_->caseOpen
            (
                rootDir_.c_str(),
                caseName_.c_str()
            );
        }
    }
    CATCH_ALL(functionName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FoamX::ICaseServerImpl::~ICaseServerImpl()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::~ICaseServerImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Delete database
    if (dbPtr_)
    {
        delete dbPtr_;
        dbPtr_ = NULL;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ICaseServerImpl::managed()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::managed()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return managed_;
}

void FoamX::ICaseServerImpl::managed(CORBA::Boolean managed)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::managed(CORBA::Boolean managed)";

    LogEntry log(functionName, __FILE__, __LINE__);

    managed_ = managed;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ICaseServerImpl::caseRoot()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::caseRoot()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return CORBA::string_dup(rootDir_.c_str());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ICaseServerImpl::caseName()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::caseName()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return CORBA::string_dup(caseName_.c_str());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IFoamProperties_ptr FoamX::ICaseServerImpl::foamProperties()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::foamProperties()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Return a reference to the TypeDescriptor object.
    IFoamProperties_ptr pProperties = IFoamProperties::_nil();

    if (!CORBA::is_nil(foamProperties_))
    {
        pProperties = IFoamProperties::_duplicate(foamProperties_);
    }

    return pProperties;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IApplication_ptr FoamX::ICaseServerImpl::application()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::application()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Return a reference to the Application object.
    IApplication_ptr pAppClass = IApplication::_nil();

    if (static_cast<FoamXServer::CaseServer::IApplication*>(app_) != NULL)
    {
        pAppClass = IApplication::_duplicate(app_);
    }

    return pAppClass;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::ICaseServerImpl::availableTimeSteps()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::availableTimeSteps()";

    LogEntry log(functionName, __FILE__, __LINE__);

    StringList* timeStepList = NULL;

    try
    {
        // Get list of time directories.
        instantList times = time().times();

        // Copy time steps into a new string list and return.
        timeStepList = new StringList();
        timeStepList->length(times.size());
        forAll(times, i)
        {
            (*timeStepList)[i] = times[i].name().c_str();
        }
    }
    CATCH_ALL(functionName);

    return timeStepList;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ICaseServerImpl::getTime()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::getTime()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Duplicate and return.
        return CORBA::string_dup(time().timeName().c_str());
    }
    CATCH_ALL(functionName);

    return CORBA::string_dup("");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::setTime
(
    const char* timeName,
    const CORBA::Long timeIndex
)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::setTime(const char*, const CORBA::Long)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        scalar timeValue(readScalar(IStringStream(timeName)()));

        time().setTime(timeValue, timeIndex);
    }
    CATCH_ALL(functionName);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ICaseServerImpl::meshDefined()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::meshDefined()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Return true if a mesh description dictionary
    fileName meshDict =
        rootDir_/caseName_/fileName("constant/polyMesh/boundary");

    return exists(meshDict);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::readMesh()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::readMesh()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // See if the mesh has been defined.
        if (!meshDefined())
        {
            throw FoamXError
            (
                E_FAIL,
                "Mesh description dictionary not found",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Create a list of entries from the boundary file.
        IOPtrList<entry> patchEntries
        (
            IOobject
            (
                "boundary",
                time().constant(),
                "polyMesh",
                time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        log << "Deleting all patches" << endl;

        // Delete all previous patch information.
        deleteAllPatches();

        log << "Adding new patches.." << endl;

        // Add patch objects for each patch defined in the mesh.
        forAll(patchEntries, i)
        {
            // Read boundary type for this patch from the
            // mesh description dictionary.
            word patchName(patchEntries[i].keyword());
            log << "Adding patch " << patchName << ".." << endl;

            // Create a new patch descriptor object for this patch
            // and load patch data from the dictionary.
            PatchProperties* patchProps = new PatchProperties(patchName);

            if (patchProps == NULL)
            {
                throw FoamXError
                (
                    E_INVALID_ARG,
                    "Failed to create PatchProperties object",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            patchProps->load(patchEntries[i].dict());

            // Add new patch to list.
            addPatch(patchProps);
        }

        // Read the field values from the field dictionaries.
        StringList_var fieldNames = app_->fields(); // Auto-release.
        for (unsigned int i = 0; i < fieldNames->length(); i++)
        {
            word fieldName(fieldNames[i]);

            // See if the field dictionary exists.
            fileName fieldDictName = time().timePath()/fieldName;

            if (exists(fieldDictName))
            {
                // Get the field object.
                IGeometricFieldImpl* pFieldValue = fieldValueMap_[fieldName];

                // Read the values from the field dictionary.
                dictionary fieldDict((IFstream(fieldDictName)()));
                pFieldValue->load(fieldDict);
            }
        }


        // Correct all the field boundary conditions for changes in mesh
        for
        (
            Dictionary<PatchProperties>::iterator iter = patchMap_.begin();
            iter != patchMap_.end();
            ++iter
        )
        {
            ///Info<< "Trying to match fields to " << iter().physicalType()
            //    << endl;
            //bool matches = fieldsMatchPatchPhysicalType(iter().physicalType());
            //Info<< "Matches:" << matches << endl;

            setPatchPhysicalType
            (
                iter().patchName().c_str(),
                iter().physicalType().c_str()
            );

            // In sync with file
            iter().modified(false);
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::importMesh
(
    const char* hostName,
    const char* rootDir,
    const char* caseName
)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::importMesh(const char* hostName, "
        "const char* rootDir, const char* caseName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName importMeshDir = fileName(rootDir)/caseName/"constant/polyMesh";
        fileName destConstDir = rootDir_/caseName_/"constant";

        if (hostName == Foam::hostName())
        {
            if (dir(importMeshDir))
            {
                cp(importMeshDir, destConstDir);
            }
            else
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Source directory " + importMeshDir
                  + " does not exist or is unreadable",
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }
        else
        {
            stringList argsList = procControl_.remoteCpArgs
            (
                userName(),
                hostName,
                importMeshDir,
                destConstDir
            );

            if (procControl_.system(argsList) != 0)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Error while copying directory " + importMeshDir + " to " 
                   + destConstDir + " using remote copy command "
                   + procControl_.commandString(argsList),
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::getFieldValues
(
    const char* fieldName,
    IGeometricField_out fieldValues
)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::getFieldValues"
        "(const char* fieldName, IGeometricField_out fieldValues)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Check for a valid field name.
        if (!fieldValueMap_.found(fieldName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid field name",
                functionName,
                __FILE__, __LINE__
            );
        }

        fieldValues = fieldValueMap_[fieldName]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::ICaseServerImpl::patchNames()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::patchNames()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return the patch name list.
    return new StringList
    (
        FoamXWordList(static_cast<const wordList&>(patchMap_.toc()))
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::addPatch(PatchProperties* patchProps)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::addPatch(PatchProperties* patchProps)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Get the patch name and the boundary type.
        // patchPhysicalType is the name of the boundary type name, not the
        // definition key.
        const word& patchName    = patchProps->patchName();
        const word& patchPhysicalType = patchProps->physicalType();

        log << "Patch = " << patchName << endl;

        // Check patch name.
        if (patchMap_.found(patchName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid patch name",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Get the specified boundary type descriptor object.
        // If the boundary type name is not valid, an exception will be thrown.
        IPatchPhysicalTypeDescriptor_var boundaryDescriptor;
        app_->findPatchPhysicalType
        (
            patchPhysicalType.c_str(),
            boundaryDescriptor.out()
        );

        // Make sure that the boundary type name is valid.
        if (CORBA::is_nil(boundaryDescriptor))
        {
            throw FoamXError
            (
                E_FAIL,
                "Invalid boundary type name '" + patchPhysicalType + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Append patch descriptor to list.
        patchMap_.append(patchName, patchProps);

        // Get patch field types from boundary descriptor.
        StringPairList_var patchFieldTypes =
            boundaryDescriptor->patchFieldTypes();

        // Add this patch to all field value objects.
        for (unsigned int i = 0; i <patchFieldTypes->length(); i++)
        {
            word fieldName(patchFieldTypes[i].name);

            if (fieldValueMap_.found(fieldName))
            {
                IGeometricFieldImpl* pFieldValue = fieldValueMap_[fieldName];

                // Convert between the patch field type key and the Foam name.

                // Get the PatchFieldDescriptor object.
                ITypeDescriptor_var pfDescriptor;
                foamProperties_->getPatchFieldType
                (
                    patchFieldTypes[i].value,
                    pfDescriptor.out()
                );

                // Make sure that the patch field type name is valid.
                if (CORBA::is_nil(pfDescriptor))
                {
                    throw FoamXError
                    (
                        E_FAIL,
                        "Invalid patch field type key '" 
                      + word(patchFieldTypes[i].value)
                      + "' for boundary type '" + patchPhysicalType + "'",
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                // Get the Foam name and set the patch field type.
                CORBA::String_var patchFieldType = pfDescriptor->name();
                pFieldValue->addPatch(patchName.c_str(), patchFieldType);
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::deletePatch(const char* patchName)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::deletePatch(const char* patchName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Check patch name.
        if (!patchMap_.found(patchName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid patch name",
                functionName,
                __FILE__, __LINE__
            );
        }
        patchMap_.remove(patchName);

        // Remove this patch from all field value objects.
        StringList_var fieldNames = app_->fields();  // Auto-release.
        for (unsigned int i=0; i <fieldNames->length(); i++)
        {
            word fieldName(fieldNames[i]);
            if (fieldValueMap_.found(fieldName))
            {
                IGeometricFieldImpl* pFieldValue = fieldValueMap_[fieldName];
                pFieldValue->deletePatch(patchName);
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::deleteAllPatches()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::deleteAllPatches()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Delete all currently defined patch information.
        wordList patchNames(patchMap_.toc());
        forAll(patchNames, i)
        {
            deletePatch(patchNames[i].c_str());
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::setPatchPhysicalType
(
    const char* patchName,
    const char* patchPhysicalType
)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::setPatchPhysicalType"
        "(const char* patchName, const char* patchPhysicalType)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Check patch name.
        if (!patchMap_.found(patchName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid patch name",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Get the specified boundary type descriptor object.
        // If the boundary type name is not valid, an exception will be thrown.
        IPatchPhysicalTypeDescriptor_var boundaryDescriptor;
        app_->getPatchPhysicalType(patchPhysicalType, boundaryDescriptor.out());

        // Make sure that the boundary type name is valid.
        if (CORBA::is_nil(boundaryDescriptor))
        {
            throw FoamXError
            (
                E_FAIL,
                "Invalid boundary type name '" + word(patchPhysicalType) + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Update the patch properties object.
        PatchProperties* patchProps = patchMap_.lookup(patchName);

        patchProps->patchType(boundaryDescriptor->patchType());
        patchProps->physicalType(patchPhysicalType);

        // Get patch field types from boundary descriptor.
        StringPairList_var patchFieldTypes =
            boundaryDescriptor->patchFieldTypes();

        // Loop over all fields and update the field value objects.
        for (unsigned int i = 0; i <patchFieldTypes->length(); i++)
        {
            word fieldName(patchFieldTypes[i].name);

            if (fieldValueMap_.found(fieldName))
            {
                // Get the patch field type descriptor.
                ITypeDescriptor_var pfDescriptor;
                foamProperties_->getPatchFieldType
                (
                    patchFieldTypes[i].value,
                    pfDescriptor.out()
                );

                // Need to pass the foam name, not the key,
                // for the patch field type.
                CORBA::String_var pfName = pfDescriptor->name();

                // Get the field object.
                IGeometricFieldImpl* pFieldValue = fieldValueMap_[fieldName];

                pFieldValue->setPatchFieldType(patchName, pfName);
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::getPatchPhysicalType
(
    const char* patchName,
    CORBA::String_out patchPhysicalType
)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::getPatchPhysicalType"
        "(const char* patchName, CORBA::String_out patchPhysicalType)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Patch = " << patchName << endl;

        // Check for a valid patch name.
        if (!patchMap_.found(patchName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid patch name",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return the boundary type for this patch.
        PatchProperties* patchProps = patchMap_.lookup(patchName);

        patchPhysicalType = CORBA::string_dup(patchProps->physicalType().c_str());
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::getDictionary
(
    const char* dictionaryName,
    CORBA::Boolean forceRead,
    IDictionaryEntry_out dictRoot
)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::getDictionary"
        "(const char*, CORBA::Boolean, IDictionaryEntry_out)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        RootDictionary* pDictRoot = NULL;

        // Force rereading just by clearing the cache.
        if (forceRead)
        {
            ObjRefHashTable<RootDictionary*>::iterator iter =
                dictionaryMap_.find(dictionaryName);

            dictionaryMap_.erase(iter);
        }

        // See if we have this dictionary entry object cached.
        if (dictionaryMap_.found(dictionaryName))
        {
            log << "Dictionary " << dictionaryName << " existing" << endl;

            pDictRoot = dictionaryMap_[dictionaryName];
        }
        else
        {
            log << "Creating new dictionary " << dictionaryName << endl;

            // Create and initialise a new dictionary entry object.
            ITypeDescriptor_var typeDescRef; // Auto release.
            app_->getDictionary(dictionaryName, typeDescRef.out());

            // Create a default RootDictionary object for this dictionary.
            pDictRoot = new RootDictionary
            (
                typeDescRef,
                rootDir_,
                caseName_
            );
            if (pDictRoot == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Couldn't create dictionary entry object",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Load the values from the file if it exists
            pDictRoot->load();

            // Add to dictionary map.
            dictionaryMap_.insert(dictionaryName, pDictRoot);
        }

        // Return a reference to the dictionary entry object.
        dictRoot = pDictRoot->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Long FoamX::ICaseServerImpl::fileModificationDate
(
    const char* dictionaryName
)
{
    fileName fName(rootDir_/caseName_/dictionaryName);

    if (exists(fName))
    {
        return label(lastModified(fName));
    }
    else
    {
        // Ugly. Let's hope not many dates with time_t==0
        return 0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::readFile
(
    const char* fName,
    CORBA::String_out contents
)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::readRawDictionary"
        "(const char* fName, CORBA::String_out contents)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName dictName = rootDir_/caseName_/fName;

        log << "File = " << dictName << endl;

        // If file exists, return its contents.
        if (exists(dictName))
        {
            IFstream ifFile(dictName);

            if (ifFile)
            {
                OStringStream fileString;

                char ch;
                while (ifFile.get(ch))
                {
                    fileString << ch;
                }

                // Must duplicate string to return.
                contents = CORBA::string_dup(fileString.str().c_str());
            }
            else
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to read file '" + dictName + "'",
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to open file '" + dictName + "'",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::writeFile
(
    const char* fName,
    const char* contents
)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::writeRawDictionary"
        "(const char* fName, const char* contents)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName dictName = rootDir_/caseName_/fName;

        log << "File = " << dictName << endl;

        // Check if dictionary already exists, and if so move to .bak
        if (file(dictName))
        {
            mv(dictName, dictName + ".bak");
        }

        // Make sure path exists.
        if (!dir(dictName.path()) && !mkDir(dictName.path()))
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to create directory '" + dictName.path() + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Write dictionary contents.
        OFstream ofFile(dictName);
        if (ofFile)
        {
            ofFile << contents;
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to write file '" + dictName + "'",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Long FoamX::ICaseServerImpl::runCase(const char* arguments)
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::runCase(const char* arguments)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Convert the arguments into a stringList and start
        // the calculation process.

        string stringArgs(arguments);

        int nArgs = 1;
        if (stringArgs.size() != 0)
        {
            nArgs += 1 + stringArgs.count(' ');
        }
        stringList args(nArgs + 1);
        args[0] = word(app_->name());

        int argi = 1;
        args[argi++] = "-case " + rootDir_/caseName_;

        size_t startPos = 0;
        int endPos;
        while ((endPos = stringArgs.find(" ", startPos)) != -1)
        {
            string oneArg(stringArgs(startPos, (endPos - startPos)));
            args[argi++] = oneArg;
            startPos = endPos + 1;
        }

        if ((stringArgs.size() != 0) && (startPos != (stringArgs.size() - 1)))
        {
            string oneArg
            (
                stringArgs
                (
                    startPos, 
                    stringArgs.size() - 1 - startPos
                )
            );
            args[argi++] = oneArg;
        }

        pid_ = procControl_.fork(args, "");

        if (pid_ == -1)
        {
            throw FoamXError
            (
                E_FAIL,
                "Error executing foam calculation process",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);

    return pid_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::killCase()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::killCase()";

    LogEntry log(functionName, __FILE__, __LINE__);

    if (pid_ > 0)
    {
        if (procControl_.kill(pid_) != 0)
        {
            throw FoamXError
            (
                E_FAIL,
                string("Error killing job with PID ") + name(pid_),
                functionName,
                __FILE__, __LINE__
            );
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::validate()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::validate()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ICaseServerImpl::modified()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::modified()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        for
        (
            ObjRefHashTable<IGeometricFieldImpl*>::iterator iter =
                fieldValueMap_.begin();
            iter != fieldValueMap_.end();
            ++iter
        )
        {
            if (iter()->modified()) return true;
        }

        for
        (
            ObjRefHashTable<RootDictionary*>::iterator iter =
                dictionaryMap_.begin();
            iter != dictionaryMap_.end();
            ++iter
        )
        {
            if (iter()->modified()) return true;
        }

        for
        (
            Dictionary<PatchProperties>::iterator iter = patchMap_.begin();
            iter != patchMap_.end();
            ++iter
        )
        {
            if (iter().modified()) return true;
        }
    }
    CATCH_ALL(functionName);

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::save()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::save()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Save (default) controlDict
        saveControlDict();

        // Delete time and reconstruct to force reading of startTime
        // setting (Time::read does not reread startTime)
        if (dbPtr_)
        {
            delete dbPtr_;
            dbPtr_ = NULL;
        }
        dbPtr_ = new Time(Time::controlDictName, rootDir_, caseName_);


        // Make sure that the time directory exists.
        if (!dir(time().timePath()))
        {
            log << "Creating time directory " << time().timePath() << endl;

            if (!mkDir(time().timePath()))
            {
                // Throw a wobbler.
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create time directory",
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }

        // For dictionaries that do not exist, make sure that a default
        // dictionary is cached.
        StringList_var dictNames = app_->dictionaries();
        for (unsigned int i = 0; i <dictNames->length(); i++)
        {
            word dictName(dictNames[i]);

            // Get dictionary type descriptor.
            ITypeDescriptor_var typeDescRef; // Auto release.
            app_->getDictionary(dictName.c_str(), typeDescRef.out());

            // Get dictionary full file name
            // (using the path information in the type descriptor).
            CORBA::String_var name = typeDescRef->name();

            fileName dictPathName = RootDictionary::path
            (
                rootDir_,
                caseName_,
                typeDescRef->dictionaryPath()
            )/fileName(name);

            if (!exists(dictPathName))
            {
                IDictionaryEntry_var dict; // Auto Release.
                getDictionary(name, false, dict.out());
            }
        }

        // Now save the case's dictionaries. 
        // Only need to save the dictionaries in the cache.
        for
        (
            HashTable<RootDictionary*>::iterator iter =
                dictionaryMap_.begin();
            iter != dictionaryMap_.end();
            ++iter
        )
        {
            // Validate the dictionary before we save.
            try
            {
                iter()->validate();
                iter()->save();
            }
            catch (ValidationError vex)
            {
                log << "Dictionary validation failed : "
                    << vex.errorMessage << " at path '" << vex.itemPath
                    << "'" << endl;
            }
        }


        // Save the patch and boundary condition information.
        writePatchData();

        // Get an ordered list of the patch names
        // Used to ensure the patchFields are written in the same order
        wordList patchNames(patchMap_.size());
        label i = 0;
        for
        (
            Dictionary<PatchProperties>::iterator iter = patchMap_.begin();
            iter != patchMap_.end();
            ++iter
        )
        {
            patchNames[i++] = iter().patchName();
        }

        // Save field dictionaries.
        StringList_var fieldNames = app_->fields();   // Auto-release.
        for (unsigned int i = 0; i <fieldNames->length(); i++)
        {
            word fieldName(fieldNames[i]);

            // Get field value object.
            IGeometricFieldImpl* pFieldValue = fieldValueMap_[fieldName];

            // Create dictionary writer object.
            DictionaryWriter dictWriter
            (
                rootDir_,
                caseName_,
                time().timeName(),
                fieldName
            );

            pFieldValue->save(dictWriter, patchNames);
        }

        log << "Case " << time().path() << " saved successfully" << endl;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseServerImpl::close()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::close()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        caseBrowser_->unlockCase(rootDir_.c_str(), caseName_.c_str());

        // Unregister object from the naming service.
        NameServer fxNameServer(orb_.orbPtr_);
        fxNameServer.unbindObject(objectName_);

        // Remove case context
        fxNameServer.removeContext
        (
            hostContext_/userContext_/string::validate<word>(rootDir_)
        );

        // Remove user context
        fxNameServer.removeContext(hostContext_/userContext_);

        fxNameServer.disconnect();

        // Delete database
        if (dbPtr_)
        {
            delete dbPtr_;
            dbPtr_ = NULL;
        }

        // Shutdown th'Orb.
        // Removed because it screws-up Java 1.4 CORBA connection
        orb_.orbPtr_->shutdown(false);       // Do not wait for completion.
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //

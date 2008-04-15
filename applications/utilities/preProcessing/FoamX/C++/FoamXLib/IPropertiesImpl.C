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

// Foam header files.
#include "IFstream.H"
#include "OSspecific.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "IPropertiesImpl.H"
#include "ITypeDescriptorImpl.H"
#include "IApplicationImpl.H"
#include "IDictionaryEntryImpl.H"
#include "RootDictionary.H"
#include "IGeometryDescriptorImpl.H"
#include "IPatchDescriptorImpl.H"
#include "ITypeDescriptorImpl.H"
#include "LogEntry.H"
#include "Paths.H"

// Namespaces
#include "FoamXNameSpaces.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::ApplicationDescriptor*
FoamX::IPropertiesImpl::readApplicationDescriptor
(
    const word& name,
    const fileName& category,
    const fileName& path,
    const bool systemClass
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::readApplicationDescriptor"
        "(const word& name, const fileName& category, "
        "const bool systemClass)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Create and initialise a new ApplicationDescriptor object.
        ApplicationDescriptor* pDesc = new ApplicationDescriptor();

        pDesc->name = name.c_str();
        pDesc->category = category.c_str();
        pDesc->path = path.c_str();
        pDesc->systemClass = systemClass;

        return pDesc;
    }
    CATCH_ALL(functionName);
}


void FoamX::IPropertiesImpl::readApplicationDescriptors
(
    Foam::HashPtrTable<FoamXServer::ApplicationDescriptor>&
        appDescriptorsMap,
    const fileName& dir,
    const fileName& category,
    const bool systemClass
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::readApplicationDescriptor"
        "(const fileName& dir, "
        "const fileName& category, "
        "const bool systemClass)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileNameList appDirs = readDir(dir, fileName::DIRECTORY);

        forAll (appDirs, i)
        {
            fileName path = dir/appDirs[i];
            fileName appDictPathName =
                path/"FoamX"/(appDirs[i] + ".cfg");

            if
            (
                file(appDictPathName)
            && !appDescriptorsMap.found(appDirs[i])
            )
            {
                appDescriptorsMap.insert
                (
                    appDirs[i],
                    readApplicationDescriptor
                    (
                        appDirs[i],
                        category,
                        path,
                        systemClass
                    )
                );
            }

            readApplicationDescriptors
            (
                appDescriptorsMap,
                path,
                category/appDirs[i],
                systemClass
            );
        }
    }
    CATCH_ALL(functionName);
}


void FoamX::IPropertiesImpl::addPatchFields
(
    const dictionary& patchFieldDict
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::addPatchFields"
        "(const dictionary& patchFieldDict)";

    // Initialise the patch field descriptor objects.
    forAllConstIter(dictionary, patchFieldDict, iter)
    {
        const word& patchFieldTypeName = iter().keyword();

        // Create and initialise a new PatchFieldDescriptor object.
        ITypeDescriptorImpl* pPatchFieldDescriptor = new ITypeDescriptorImpl
        (
            patchFieldTypeName, 
            patchFieldDict.name(),
            iter(),
            foamTypesDict_
        );

        if (pPatchFieldDescriptor == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to create PatchFieldDescriptor object for "
                "patch field type '" + patchFieldTypeName + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Add to map.
        patchFieldMap_.insert(patchFieldTypeName, pPatchFieldDescriptor);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IPropertiesImpl::IPropertiesImpl(bool readOnly)
:
    readOnly_(readOnly)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::IPropertiesImpl(bool readOnly)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Open the FoamX system config dictionary.
        fileName fxConfigFileName = dotFoam("apps/FoamX/FoamX.cfg");
        if (!exists(fxConfigFileName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot find FoamX root configuration dictionary '"
              + fxConfigFileName + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }
        dictionary configDict((IFstream(fxConfigFileName)()));

        // Read list data from config dictionary.
        availableModules_.read(configDict.lookup("availableModules"));


        // Open the Foam types dictionary.
        fileName foamTypesDictFileName = Paths::config/"types/types.cfg";

        if (!exists(foamTypesDictFileName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot find Foam types dictionary '"
              + foamTypesDictFileName + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }
        (IFstream(foamTypesDictFileName)()) >> foamTypesDict_;

        // Initialise the foam type descriptor objects.
        forAllConstIter(dictionary, foamTypesDict_, iter)
        {
            const word& foamTypeName = iter().keyword();

            // Create and initialise a new TypeDescriptor object.
            ITypeDescriptorImpl* pFieldTypeDescriptor = new ITypeDescriptorImpl
            (
                foamTypeName,
                "FoamX:foamTypes",
                iter(),
                foamTypesDict_
            );

            if (pFieldTypeDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create Foam TypeDescriptor object for "
                    "field type '" + foamTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Add to map.
            foamTypeMap_.insert(foamTypeName, pFieldTypeDescriptor);
        }


        dictionary geometryDict
        (
            IFstream(Paths::config/"types/geometries.cfg")()
        );

        // Initialise the geometry descriptor objects.
        forAllConstIter(dictionary, geometryDict, iter)
        {
            const word& geometricTypeName = iter().keyword();

            // Create and initialise a new GeometryDescriptor object.
            IGeometryDescriptorImpl* pGeometryTypeDescriptor = new
            IGeometryDescriptorImpl
            (
                geometricTypeName,
                geometryDict
            );

            if (pGeometryTypeDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create field GeometryDescriptor object for "
                    "geometry type '" + geometricTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Add to map.
            geometryTypeMap_.insert(geometricTypeName, pGeometryTypeDescriptor);
        }

        dictionary patchDict(IFstream(Paths::config/"types/patches.cfg")());

        // Initialise the patch descriptor objects.
        forAllConstIter(dictionary, patchDict, iter)
        {
            const word& patchTypeName = iter().keyword();

            // Create and initialise a new PatchDescriptor object.
            IPatchDescriptorImpl* pPatchDescriptor =
                new IPatchDescriptorImpl(patchTypeName, patchDict);

            if (pPatchDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create PatchDescriptor object for patch type '"
                   + patchTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Add to map.
            patchMap_.insert(patchTypeName, pPatchDescriptor);
        }


        dictionary patchFieldDict
        (
            IFstream(Paths::config/"types/patchFields.cfg")()
        );

        addPatchFields(patchFieldDict);


        // ---------------------------------------------------------------------
        // Open the user's Foam control dictionary.
        fileName controlDictFileName = dotFoam("controlDict");
        if (exists(controlDictFileName))
        {
            dictionary controlDict((IFstream(controlDictFileName)()));

            stringList rawRootDirs(1, ".");

            if (controlDict.found("caseRoots"))
            {
                // Read root directory list.
                // Filter out directories which do not exist.
                stringList controlDictRootDirs(controlDict.lookup("caseRoots"));

                if (controlDictRootDirs.size())
                {
                    rawRootDirs = controlDictRootDirs;
                }
            }

            // Check that each root directory actually exists.
            forAll (rawRootDirs, i)
            {
                fileName rootDir = rawRootDirs[i];
                rootDir.expand();

                if (exists(rootDir))
                {
                    rootDirectories_.length(rootDirectories_.length() + 1);
                    rootDirectories_[rootDirectories_.length() - 1]
                        = rootDir.c_str();

                    rawRootDirectories_.length
                    (
                        rawRootDirectories_.length() + 1
                    );

                    rawRootDirectories_[rawRootDirectories_.length() - 1]
                        = rawRootDirs[i].c_str();
                }
                else
                {
                    WarningIn(functionName)
                        << "User specified root directory " << rawRootDirs[i];

                    if (rootDir.size() != rawRootDirs[i].size())
                    {
                        Info<< " (" << rootDir << ")";
                    }

                    Info<< " not found. Removing from root directories list."
                        << endl;
                }
            }
        }


        // ---------------------------------------------------------------------
        // Initialise the application class descriptor map.

        // Scan the user's utilities directory for valid solvers
        readApplicationDescriptors
        (
            appDescriptorMap_,
            Paths::userSolvers,
            "",
            false
        );

        // Scan the system utilities directory for valid solvers
        readApplicationDescriptors
        (
            appDescriptorMap_,
            Paths::solvers,
            "",
            true
        );


        // ---------------------------------------------------------------------
        // Initialise the utility descriptor map.

        // Scan the user's utilities directory for valid utilities
        readApplicationDescriptors
        (
            utilityDescriptorMap_,
            Paths::userUtilities,
            "",
            false
        );

        // Scan the system utilities directory for valid utilities
        readApplicationDescriptors
        (
            utilityDescriptorMap_,
            Paths::utilities,
            "",
            true
        );
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IPropertiesImpl::~IPropertiesImpl()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::~IPropertiesImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::availableModules()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::availableModules()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(availableModules_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::rootDirectories()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::rootDirectories()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(rootDirectories_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::rawRootDirectories()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::rawRootDirectories()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(rawRootDirectories_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::addRootDirectory(const char* rawRootDir)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::addRootDirectory(const char* rawRootDir)";

    LogEntry log(functionName, __FILE__, __LINE__);

    fileName rootDir = rawRootDir;
    rootDir.expand();    

    if (!rootDirectories_.found(rootDir.c_str()))
    {
        rootDirectories_.append(rootDir.c_str());
        rawRootDirectories_.append(rawRootDir);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::deleteRootDirectory(const char* rawRootDir)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::deleteRootDirectory(const char* rawRootDir)";

    LogEntry log(functionName, __FILE__, __LINE__);

    fileName rootDir = rawRootDir;
    rootDir.expand();    

    if (rootDirectories_.found(rootDir.c_str()))
    {
        rootDirectories_.remove(rootDir.c_str());
        rawRootDirectories_.remove(rawRootDir);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::foamTypes()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::foamTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return
    return new FoamXServer::StringList
    (
        FoamXWordList(static_cast<const wordList&>(foamTypeMap_.toc()))
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::getFoamType
(
    const char* foamTypeKey,
    ITypeDescriptor_out typeDesc
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getFoamType"
        "(const char* foamTypeKey, ITypeDescriptor_out typeDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // Check for valid foam type name.
        if (!foamTypeMap_.found(foamTypeKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid foam type name '" + word(foamTypeKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to the TypeDescriptor object.
        typeDesc = foamTypeMap_[foamTypeKey]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::geometryTypes()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::geometryTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return
    return new FoamXServer::StringList
    (
        FoamXWordList(static_cast<const wordList&>(geometryTypeMap_.toc()))
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::getGeometryType
(
    const char* geometryTypeKey,
    IGeometryDescriptor_out geometryDesc
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getGeometryType"
        "(const char* geometryTypeKey, IGeometryDescriptor_out geometryDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // Check for valid geometry type name.
        if (!geometryTypeMap_.found(geometryTypeKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid geometry type name '" + word(geometryTypeKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to the GeometryDescriptor object.
        geometryDesc = geometryTypeMap_[geometryTypeKey]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::patchTypes()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::patchTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return
    return new FoamXServer::StringList
    (
        FoamXWordList(static_cast<const wordList&>(patchMap_.toc()))
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void  FoamX::IPropertiesImpl::getPatchType
(
    const char* patchTypeKey,
    IPatchDescriptor_out patchDesc
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getPatchType"
        "(const char* patchTypeKey, IPatchDescriptor_out patchDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // Check for valid patch name.
        if (!patchMap_.found(patchTypeKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid patch type name '" + word(patchTypeKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to the PatchDescriptor object.
        patchDesc = patchMap_[patchTypeKey]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void  FoamX::IPropertiesImpl::findPatchType
(
    const char* patchTypeName,
    IPatchDescriptor_out patchDesc
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::findPatchType"
        "(const char* patchTypeName, IPatchDescriptor_out patchDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        IPatchDescriptorImpl* pPatchDescriptor = NULL;

        if (patchMap_.found(patchTypeName));
        {
            pPatchDescriptor = patchMap_.find(patchTypeName)();
        }

        // Return a reference to the PatchFieldDescriptor object.
        if (pPatchDescriptor != NULL)
        {
            patchDesc = pPatchDescriptor->_this();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IPropertiesImpl::patchFieldTypes()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::patchFieldTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return
    return new FoamXServer::StringList
    (
        FoamXWordList(static_cast<const wordList&>(patchFieldMap_.toc()))
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::getPatchFieldType
(
    const char* patchFieldTypeKey,
    ITypeDescriptor_out patchFieldDesc
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getPatchFieldType"
        "(const char* patchFieldTypeKey,"
        "ITypeDescriptor_out patchFieldDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Check for valid patch field name.
        if (!patchFieldMap_.found(patchFieldTypeKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid patch field type name '"
              + word(patchFieldTypeKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to the PatchFieldDescriptor object.
        patchFieldDesc = patchFieldMap_[patchFieldTypeKey]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::findPatchFieldType
(
    const char* patchFieldTypeName,
    ITypeDescriptor_out patchFieldDesc
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::findPatchFieldType"
        "(const char* patchFieldTypeName, "
        "ITypeDescriptor_out patchFieldDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        ITypeDescriptorImpl* pPatchFieldDescriptor = NULL;

        if (patchFieldMap_.found(patchFieldTypeName))
        {
            pPatchFieldDescriptor = patchFieldMap_.find(patchFieldTypeName)();
        }

        // Return a reference to the PatchFieldDescriptor object.
        if (pPatchFieldDescriptor != NULL)
        {
            patchFieldDesc = pPatchFieldDescriptor->_this();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::getFoamControlDict
(
    FoamXServer::IDictionaryEntry_out controlDict
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getFoamControlDict"
        "(FoamXServer::ApplicationDescriptor_out utilityDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName controlDictCfgPath = 
            Paths::config/"dictionaries/OpenFOAMControlDict/controlDict.cfg";

        dictionary controlDictConfigDict((IFstream(controlDictCfgPath)()));

        ITypeDescriptorImpl* controlDictDesc = new ITypeDescriptorImpl
        (
            "controlDict",
            controlDictCfgPath,
            controlDictConfigDict,
            foamTypesDict_
        );

        fileName controlDictPath = Foam::dotFoam("");

        // Create an appropriate sub entry object and store its reference.
        RootDictionary* controlDictPtr = new RootDictionary
        (
            controlDictDesc->_this(),
            controlDictPath,
            ""
        );

        if (controlDictPtr == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Couldn't create IDictionaryEntryImpl object for OpenFOAM "
                "controlDict " + controlDictPath,
                functionName,
                __FILE__, __LINE__
            );
        }

        controlDict = controlDictPtr->_this();

        // Load the values from the file.
        controlDictPtr->load();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::ApplicationDescriptorList* FoamX::IPropertiesImpl::applicationes()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::applicationes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Construct an application class list and return.
    ApplicationDescriptorList* pAppClassList =
        new ApplicationDescriptorList();
    pAppClassList->length(appDescriptorMap_.size());

    label i = 0;

    for
    (
        Foam::HashPtrTable<ApplicationDescriptor>::iterator iter = 
            appDescriptorMap_.begin();
        iter != appDescriptorMap_.end();
        ++iter
    )
    {
        (*pAppClassList)[i++] = *iter();
    }

    return pAppClassList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::getApplication
(
    const char* appKey,
    IApplication_out app
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getApplication"
        "(const char* appKey, IApplication_out app)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        IApplicationImpl* pAppClass = NULL;

        // See if the specified application class name is valid.
        if (!appDescriptorMap_.found(appKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "getApplication::Invalid application class name '"
              + word(appKey) + "'.",
                "IPropertiesImpl",
                __FILE__, __LINE__
            );
        }

        // See if we have this application class object cached.
        if (appMap_.found(appKey))
        {
            log << "Existing." << endl;
            pAppClass = appMap_[appKey];
        }
        else
        {
            log << "New application class " << app << endl;

            // Create and initialise a new Application object.
            pAppClass = new IApplicationImpl
            (
                *appDescriptorMap_[appKey],
                *this
            );

            if (pAppClass == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create Application object for '"
                  + word(appKey) + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Load the application class from the definition dictionary.
            // Allow user defined application classes to override the system
            // classes.
            pAppClass->load();

            // Add application class object to map.
            appMap_.insert(appKey, pAppClass);
        }

        // Return a reference to the Application object.
        app = pAppClass->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::addApplication
(
    const char* appKey,
    IApplication_out app
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::addApplication"
        "(const char* appKey, IApplication_out app)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // If this object is read-only, throw a dicky fit.
        if (readOnly_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid call to addApplication for '" + word(appKey)
              + "'. Object is read only.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if the specified application class name is valid.
        if (appDescriptorMap_.found(appKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid application class name '" + word(appKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Create and initialise a new ApplicationDescriptor object.
        ApplicationDescriptor* pDesc = new ApplicationDescriptor();
        pDesc->name        = appKey;
        pDesc->category    = "";
        pDesc->systemClass = false;

        // Create and initialise a new Application object.
        // User defined application class.
        IApplicationImpl* pAppClass = new IApplicationImpl
        (
            *pDesc,
            *this
        );

        if (pAppClass == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to create Application object for '"
              + word(appKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Add application class object to map.
        appMap_.insert(appKey, pAppClass);

        // Add application class descriptor to map.
        appDescriptorMap_.insert(appKey, pDesc);

        // Return a reference to the application class object.
        app = pAppClass->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::deleteApplication(const char* appKey)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::deleteApplication"
        "(const char* appKey)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // If this object is read-only, throw a wobbler.
        if (readOnly_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid call to deleteApplication for " 
              + word(appKey) + "'. Object is read only.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if the specified application class name is valid.
        if (!appDescriptorMap_.found(appKey))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid application class name '" + word(appKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if the specified application class name is user defined.
        if (appDescriptorMap_[appKey]->systemClass)
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Unable to delete system application class '"
              + word(appKey) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Remove the complete application directory structure
        if (!rmDir(fileName(appDescriptorMap_[appKey]->path)))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Unable to delete application class "
              + fileName(appDescriptorMap_[appKey]->path),
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if we have this application class object cached.
        if (appMap_.found(appKey))
        {
            ObjRefHashTable<IApplicationImpl*>::iterator iter
            (
                appMap_.find(appKey)
            );

            appMap_.erase(iter);   // Releases implementation reference.
        }

        // Remove application class descriptor.
        Foam::HashPtrTable<ApplicationDescriptor>::iterator iter
        (
            appDescriptorMap_.find(appKey)
        );
        appDescriptorMap_.erase(iter);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::cloneApplication
(
    const char* appKeySrc,
    const char* appKeyDest,
    const char* appDestPath,
    IApplication_out app
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::cloneApplication"
        "(const char* appKeySrc, const char* appKeyDest, "
        "IApplication_out app)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // If this object is read-only, throw a dicky fit.
        if (readOnly_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid call to cloneApplication for '"
              + word(appKeySrc) + "'. Object is read only.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if the specified source application class name is valid.
        if (!appDescriptorMap_.found(appKeySrc))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid source application class name '"
              + word(appKeySrc) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if the specified destination application class name is valid.
        if (appDescriptorMap_.found(appKeyDest))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid destination application class name '"
              + word(appKeyDest) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Copy application class definition files.
        ApplicationDescriptor* pDescriptor =
            appDescriptorMap_[appKeySrc];

        // Make application class directory.
        if (!cp(fileName(pDescriptor->path), appDestPath))
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Failed to copy directory " + fileName(pDescriptor->path)
              + " to " + appDestPath,
                functionName,
                __FILE__, __LINE__
            );
        }

        // Create and initialise a new ApplicationDescriptor object.
        // User defined application class.
        ApplicationDescriptor* pDesc = new ApplicationDescriptor();
        pDesc->name        = appKeyDest;
        pDesc->category    = pDescriptor->category;
        pDesc->path        = appDestPath;
        pDesc->systemClass = false;

        // Create and initialise a new Application object.
        IApplicationImpl* pAppClass = new IApplicationImpl
        (
            *pDesc,
            *this
        );

        if (pAppClass == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to create Application object '"
              + word(appKeyDest) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Load the application class details from the definition dictionary.
        // Allow user defined application classes to override the system
        // classes.
        pAppClass->load();

        // Add to application class map.
        appMap_.insert(appKeyDest, pAppClass);

        // Add application class descriptor to map.
        appDescriptorMap_.insert(appKeyDest, pDesc);

        // Return a reference to the application class object.
        app = pAppClass->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::ApplicationDescriptorList* FoamX::IPropertiesImpl::utilities()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::utilities()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Construct an application class list and return.
    ApplicationDescriptorList* pUtilityList = 
        new ApplicationDescriptorList();
    pUtilityList->length(utilityDescriptorMap_.size());

    label i = 0;

    for
    (
        Foam::HashPtrTable<ApplicationDescriptor>::iterator iter = 
            utilityDescriptorMap_.begin();
        iter != utilityDescriptorMap_.end();
        ++iter
    )
    {
        (*pUtilityList)[i++] = *iter();
    }

    return pUtilityList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::getUtilityControlDict
(
    const char* utilityName,
    const char* rootDir,
    const char* caseName,
    FoamXServer::IDictionaryEntry_out controlDict
)
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::getUtilityControlDict"
        "(const char* utilityName, const char* rootDir, const char* caseName,"
        "FoamXServer::ApplicationDescriptor_out utilityDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Check for valid utility name.
        if (!utilityDescriptorMap_.found(utilityName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid foam utility name '" + word(utilityName) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        ApplicationDescriptor* utilityDescriptor = 
            utilityDescriptorMap_[utilityName];

        word name(utilityDescriptor->name);
        fileName path(utilityDescriptor->path);

        fileName utilityCfgPath = path/"FoamX"/(name + ".cfg");
        dictionary utilityConfigDict((IFstream(utilityCfgPath)()));
        word controlDictName = name + "Dict";

        if (utilityConfigDict.found(controlDictName))
        {
            ITypeDescriptorImpl* controlDictDesc = new ITypeDescriptorImpl
            (
                controlDictName,
                utilityCfgPath,
                utilityConfigDict.subDict(controlDictName),
                foamTypesDict_
            );

            // Create an appropriate sub entry object and store its reference.
            RootDictionary* controlDictPtr = new RootDictionary
            (
                controlDictDesc->_this(),
                rootDir,
                caseName
            );

            if (controlDictPtr == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Couldn't create IDictionaryEntryImpl object for utility "
                  + word(utilityName),
                    functionName,
                    __FILE__, __LINE__
                );
            }

            controlDict = controlDictPtr->_this();

            // Load the values from the file.
            controlDictPtr->load();
        }
        else
        {
            controlDict = FoamXServer::IDictionaryEntry::_nil();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::validate()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::validate()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
    }
    CATCH_ALL(functionName);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::saveSystemProperties()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::saveSystemProperties()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // If this object is read-only, throw a dicky fit.
        if (readOnly_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid call to save. IPropertiesImpl object is read only",
                functionName,
                __FILE__, __LINE__
            );
        }

        {
        fileName configDictFileName = dotFoam("apps/FoamX/FoamX.cfg");
        DictionaryWriter dict(configDictFileName);

        dict.writeHeader
        (
            "FoamX System Properties.",
            "dictionary"
        );

        dict.writeEntry("availableModules", availableModules_);

        dict.startSubDict("processControl");
        dict.writeEntry("remoteShell", string("rsh"));
        dict.endSubDict();
        dict.writeBar();

        // Geometry type definitions.
        dict.writeSectionHeader("Geometry type definitions.");
        dict.startSubDict("geometryTypes");
        for
        (
            ObjRefHashTable<IGeometryDescriptorImpl*>::iterator iter =
                geometryTypeMap_.begin();
            iter != geometryTypeMap_.end();
            ++iter
        )
        {
            iter()->save(dict);
            dict.writeEndl();
        }
        dict.endSubDict();
        dict.writeBar();

        // Patch type definitions.
        dict.writeSectionHeader("Patch type definitions.");
        dict.startSubDict("patchTypes");
        for
        (
            ObjRefHashTable<IPatchDescriptorImpl*>::iterator iter =
                patchMap_.begin();
            iter != patchMap_.end();
            ++iter
        )
        {
            iter()->save(dict);
            dict.writeEndl();
        }
        dict.endSubDict();
        dict.writeBar();

        // Patch field type definitions.
        dict.writeSectionHeader("Patch field type definitions.");
        dict.startSubDict("patchFieldTypes");
        for
        (
            ObjRefHashTable<ITypeDescriptorImpl*>::iterator iter =
                patchFieldMap_.begin();
            iter != patchFieldMap_.end();
            ++iter
        )
        {
            iter()->save(dict, true);
            dict.writeEndl();
        }
        dict.endSubDict();
        dict.writeEndl();
        dict.writeEndBar();
        }



        {
        DictionaryWriter dict(Paths::config/"types/types.cfg");

        dict.writeHeader
        (
            "Primitive types.",
            "dictionary"
        );

        for
        (
            ObjRefHashTable<ITypeDescriptorImpl*>::iterator iter =
                foamTypeMap_.begin();
            iter != foamTypeMap_.end();
            ++iter
        )
        {
            iter()->save(dict, false);
            dict.writeEndl();
        }

        dict.writeEndl();
        dict.writeEndBar();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPropertiesImpl::saveUserProperties()
{
    static const char* functionName =
        "FoamX::IPropertiesImpl::saveUserProperties()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // If this object is read-only, throw a dicky fit.
        if (readOnly_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid call to save. IPropertiesImpl object is read only",
                functionName,
                __FILE__, __LINE__
            );
        }

        fileName controlDictFileName = dotFoam("controlDict");
        DictionaryWriter dict(controlDictFileName);

        dict.writeHeader
        (
            "FoamX User Properties.",
            "dictionary"
        );

        dict.writeEndl();
        dict.writeEntry("caseRoots", rootDirectories_);
        dict.writeEndl();
        dict.writeEndBar();
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //

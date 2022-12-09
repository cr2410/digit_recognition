========================================================================
    CONSOLE APPLICATION : digit_recognition Project Overview
========================================================================

AppWizard has created this digit_recognition application for you.

This file contains a summary of what you will find in each of the files that
make up your digit_recognition application.


digit_recognition.vcxproj
    This is the main project file for VC++ projects generated using an Application Wizard.
    It contains information about the version of Visual C++ that generated the file, and
    information about the platforms, configurations, and project features selected with the
    Application Wizard.

digit_recognition.vcxproj.filters
    This is the filters file for VC++ projects generated using an Application Wizard. 
    It contains information about the association between the files in your project 
    and the filters. This association is used in the IDE to show grouping of files with
    similar extensions under a specific node (for e.g. ".cpp" files are associated with the
    "Source Files" filter).

digit_recognition.cpp
    Training -
	  Program is initally trained for 25 utterance of each digit
	   observation sequence is generated and model is created for each utterance of digit.
	   once the 25 model of all utterances is created we take the average of all models
	   and obtain  new model , this obtained new model is now made inital model and process
	   is repeated again . The new average model is stored as final model for that digit.
	   Same process is repeated for each digts. 

   Testing-
      Program is tested for 5 utterance of each digit
	      observation sequence is generated and with each model 0f 1-10 probability
		  (O/lambda) is calculated using solution of problem 1 , whichever model 
		  gives the best probability corresponsding digit of that model is recognized
		  digit.

/////////////////////////////////////////////////////////////////////////////
Other standard files:

StdAfx.h, StdAfx.cpp
    These files are used to build a precompiled header (PCH) file
    named digit_recognition.pch and a precompiled types file named StdAfx.obj.

/////////////////////////////////////////////////////////////////////////////
Other notes:

AppWizard uses "TODO:" comments to indicate parts of the source code you
should add to or customize.

/////////////////////////////////////////////////////////////////////////////

# Configuration file for building Dynare Windows Installer
# Uses "NullSoft Scriptable Installer System", aka NSIS (see http://nsis.sourceforge.net)
# NSIS can run from both Windows and Linux (see "nsis" package in Debian)

# How to build the installer:
# - build the preprocessor, the MEX binaries (for Matlab 7.4, 7.5 and for Octave), and the documentation (PDF files + HTML manual)
# - run "makensis dynare.nsi" to create the installer

!define VERSION "4.svn20080629"

Name "Dynare ${VERSION}"

OutFile "dynare-${VERSION}-win32.exe"

InstallDir "c:\dynare\${VERSION}"

# Use the Modern User Interface (version 2)
!include MUI2.nsh

!define MUI_WELCOMEPAGE_TEXT "This wizard will guide you through the installation of Dynare ${VERSION}.$\n$\nDynare is distributed under the GNU General Public License (GPL) version 3.$\n$\nIf you accept the license, click Next button to continue the installation."
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!define MUI_FINISHPAGE_TEXT "Dynare is now installed.$\n$\nTo use it, you need to run Matlab (R) or Octave, and then to add directory '$INSTDIR\matlab' to the search path (with 'addpath' command).$\n$\nTo uninstall Dynare, go to the 'Add/Remove programs' menu of the Control Panel."
!define MUI_FINISHPAGE_LINK_LOCATION http://www.dynare.org
!define MUI_FINISHPAGE_LINK "Go to Dynare homepage"
!define MUI_FINISHPAGE_SHOWREADME $INSTDIR\doc
!define MUI_FINISHPAGE_SHOWREADME_TEXT "Open documentation directory"
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_WELCOME
!insertmacro MUI_UNPAGE_INSTFILES
!insertmacro MUI_UNPAGE_FINISH

!insertmacro MUI_LANGUAGE "English"

!define REGLOC "Software\Microsoft\Windows\CurrentVersion\Uninstall\Dynare ${VERSION}"
!define SMLOC "$SMPROGRAMS\Dynare ${VERSION}"

!macro DETERMINE_CONTEXT
 # Determine if we are admin or not
 # This will change the start menu directory and the registry root key (HKLM or HKLU)
 UserInfo::getAccountType
 Pop $0
 StrCmp $0 "Admin" +3
 SetShellVarContext current
 Goto +2
 SetShellVarContext all
!macroend

Section
!insertmacro DETERMINE_CONTEXT
 SetOutPath $INSTDIR
 File LICENSE.txt

 SetOutPath $INSTDIR\matlab
 File /r matlab\*.m
 File matlab\dynare_m.exe

 SetOutPath $INSTDIR\mex
 File /r mex\*.mexw32 mex\*.dll

 SetOutPath $INSTDIR\mex\octave
 File mex\octave\rcond.m mex\octave\*.mex

 SetOutPath $INSTDIR\doc
 File doc\guide.pdf doc\userguide\Dynare_UserGuide_WebBeta.pdf doc\bvar-a-la-sims.pdf

 SetOutPath $INSTDIR\doc\manual
 File doc\manual\*.html

 WriteUninstaller $INSTDIR\uninstall.exe

 # Create start menu entries
 CreateDirectory "${SMLOC}"
 CreateShortcut "${SMLOC}\Documentation.lnk" "$INSTDIR\doc"
 CreateShortcut "${SMLOC}\Uninstall.lnk" "$INSTDIR\uninstall.exe"

 # Create entry in "Add/Remove programs"
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "DisplayName" "Dynare ${VERSION}"
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "DisplayVersion" "${VERSION}"
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "InstallLocation" $INSTDIR
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "UninstallString" "$INSTDIR\uninstall.exe"
 WriteRegDWORD SHELL_CONTEXT "${REGLOC}" "NoModify" 1
 WriteRegDWORD SHELL_CONTEXT "${REGLOC}" "NoRepair" 1
SectionEnd

Section "Uninstall"
!insertmacro DETERMINE_CONTEXT

 # First delete the uninstaller
 Delete $INSTDIR\uninstall.exe
 Delete $INSTDIR\LICENSE.txt
 Rmdir /r $INSTDIR\matlab
 Rmdir /r $INSTDIR\mex
 Rmdir /r $INSTDIR\doc
 # We don't force deletion of installation directory (with /r), to avoid deleting important files
 Rmdir $INSTDIR

 # Delete start menu entries
 Rmdir /r "${SMLOC}"

 # Delete entry in "Add/Remove programs"
 DeleteRegKey SHELL_CONTEXT "${REGLOC}"
SectionEnd

# Configuration file for building Dynare Windows Installer
# Uses "NullSoft Scriptable Installer System", aka NSIS (see http://nsis.sourceforge.net)
# NSIS can be run from both Windows and Linux (see "nsis" package in Debian)

# NOTE: if you want to build from Debian, you'll need to replace /usr/share/nsis/Plugins/System.dll by the System.dll included in the windows distribution of NSIS (see http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=319999)

# How to build the installer:
# - build the preprocessor, the MEX binaries (for Matlab 7.4, 7.5 and for Octave), and the documentation (PDF files + HTML manual)
# - note that MEX files in the mex/2007a directory must have the .dll extension (and not the .mexw32 extension, which is only supported since Matlab 7.1)
# - run "makensis dynare.nsi" to create the installer
# - if there is no failure, this will create a file "dynare-VERSION-win32.exe" in the current directory

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
!define MUI_FINISHPAGE_LINK_LOCATION http://www.dynare.org
!define MUI_FINISHPAGE_LINK "Go to Dynare homepage"
!define MUI_FINISHPAGE_SHOWREADME $INSTDIR\README.txt
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
 File README.txt gpl.txt fdl.txt

 SetOutPath $INSTDIR\matlab
 File /r ..\matlab\*.m
 File ..\matlab\dynare_m.exe

 SetOutPath $INSTDIR\mex\2007a
 File ..\mex\2007a\*.dll

 SetOutPath $INSTDIR\mex\2007b
 File ..\mex\2007b\*.mexw32

 SetOutPath $INSTDIR\mex\octave
 File ..\mex\octave\rcond.m ..\mex\octave\*.mex

 SetOutPath $INSTDIR\doc
 File ..\doc\manual.pdf ..\doc\guide.pdf ..\doc\userguide\Dynare_UserGuide_WebBeta.pdf ..\doc\bvar-a-la-sims.pdf ..\doc\macroprocessor\macroprocessor.pdf ..\doc\preprocessor\preprocessor.pdf

 SetOutPath $INSTDIR\doc\manual-html
 File ..\doc\manual-html\*.html

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
 Delete $INSTDIR\README.txt
 Delete $INSTDIR\gpl.txt
 Delete $INSTDIR\fdl.txt
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

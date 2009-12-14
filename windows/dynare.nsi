# Configuration file for building Dynare Windows Installer
# Uses "NullSoft Scriptable Installer System", aka NSIS (see http://nsis.sourceforge.net)
# NSIS can be run from both Windows and Linux (see "nsis" package in Debian)

# NOTE: if you want to build from Debian, you'll need to replace /usr/share/nsis/Plugins/System.dll by the System.dll included in the windows distribution of NSIS (see http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=319999)

# How to build the installer:
# - build: the preprocessor, the MEX binaries (for MATLAB and for Octave), and the documentation (PDF files + HTML manual)
# - run "makensis dynare.nsi" to create the installer
# - if there is no failure, this will create a file "dynare-VERSION-win.exe" in the current directory

!include dynare-version.nsi

Name "Dynare ${VERSION}"

OutFile "dynare-${VERSION}-win.exe"

InstallDir "c:\dynare\${VERSION}"

# Use the Modern User Interface (version 2)
!include MUI2.nsh

!define MUI_WELCOMEPAGE_TEXT "This wizard will guide you through the installation of Dynare ${VERSION}.$\n$\nDynare is distributed under the GNU General Public License (GPL) version 3.$\n$\nIf you accept the license, click Next button to continue the installation."
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_COMPONENTS
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

# Strip preprocessor, dynare++ and 32-bit DLL binaries
# (64-bit DLL are compiled with MSVC and therefore are small)
!system 'strip ..\matlab\dynare_m.exe'
!system 'strip ..\dynare++\src\dynare++.exe'
!system 'strip ..\mex\matlab\win32-6.5-7.4\*'
!system 'strip ..\mex\matlab\win32-7.5-7.9\*'

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

Section "Dynare core (preprocessor and M-files)"
 SectionIn RO
!insertmacro DETERMINE_CONTEXT
 SetOutPath $INSTDIR
 File README.txt mexopts.bat ..\license.txt

 SetOutPath $INSTDIR\matlab
 File /r ..\matlab\*.m
 File ..\matlab\dynare_m.exe

 WriteUninstaller $INSTDIR\uninstall.exe

 # Create start menu entries
 CreateDirectory "${SMLOC}"
 CreateShortcut "${SMLOC}\Uninstall.lnk" "$INSTDIR\uninstall.exe"

 # Create entry in "Add/Remove programs"
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "DisplayName" "Dynare ${VERSION}"
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "DisplayVersion" "${VERSION}"
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "InstallLocation" $INSTDIR
 WriteRegStr SHELL_CONTEXT "${REGLOC}" "UninstallString" "$INSTDIR\uninstall.exe"
 WriteRegDWORD SHELL_CONTEXT "${REGLOC}" "NoModify" 1
 WriteRegDWORD SHELL_CONTEXT "${REGLOC}" "NoRepair" 1
SectionEnd

SectionGroup "MEX files for MATLAB"

Section "MEX files for MATLAB 32-bit, version 6.5 to 7.4 (R13 to R2007a)"
 SetOutPath $INSTDIR\mex\matlab\win32-6.5-7.4
 File ..\mex\matlab\win32-6.5-7.4\*.dll
SectionEnd

Section "MEX files for MATLAB 32-bit, version 7.5 to 7.9 (R2007b to R2009b)"
 SetOutPath $INSTDIR\mex\matlab\win32-7.5-7.9
 File ..\mex\matlab\win32-7.5-7.9\*.mexw32
SectionEnd

Section "MEX files for MATLAB 64-bit, version 7.2 (R2006a)"
 SetOutPath $INSTDIR\mex\matlab\win64-7.2
 File ..\mex\matlab\win64-7.2\*.mexw64
SectionEnd

Section "MEX files for MATLAB 64-bit, version 7.3 to 7.4 (R2006b to R2007a)"
 SetOutPath $INSTDIR\mex\matlab\win64-7.3-7.4
 File ..\mex\matlab\win64-7.3-7.4\*.mexw64
SectionEnd

Section "MEX files for MATLAB 64-bit, version 7.5 to 7.7 (R2007b to R2008b)"
 SetOutPath $INSTDIR\mex\matlab\win64-7.5-7.7
 File ..\mex\matlab\win64-7.5-7.7\*.mexw64
SectionEnd

Section "MEX files for MATLAB 64-bit, version 7.8 to 7.9 (R2009a to R2009b)"
 SetOutPath $INSTDIR\mex\matlab\win64-7.8-7.9
 File ..\mex\matlab\win64-7.8-7.9\*.mexw64
SectionEnd

SectionGroupEnd

Section "MEX files for Octave 3.2.3 (MinGW build)"
 SetOutPath $INSTDIR\mex\octave
 File ..\mex\octave\*.mex
SectionEnd

Section "Dynare++ (standalone executable)"
 SetOutPath $INSTDIR\dynare++
 File ..\dynare++\src\dynare++.exe
 # The list of DLLs given here is used when building with Octave/MinGW build 3.2.3 (disabling POSIX threads)
 File ..\dynare++\src\atlas.dll ..\dynare++\src\blas.dll ..\dynare++\src\cblas.dll ..\dynare++\src\lapack.dll ..\dynare++\src\libgcc_s_dw2-1.dll ..\dynare++\src\libgfortran-3.dll
SectionEnd

Section "Documentation (Dynare and Dynare++)"
 SetOutPath $INSTDIR\doc
 File ..\doc\manual.pdf ..\doc\guide.pdf ..\doc\userguide\UserGuide.pdf ..\doc\bvar-a-la-sims.pdf ..\doc\macroprocessor\macroprocessor.pdf ..\doc\preprocessor\preprocessor.pdf

 SetOutPath $INSTDIR\doc\manual-html
 File ..\doc\manual-html\*.html

 SetOutPath $INSTDIR\doc\dynare++
 File ..\dynare++\doc\dynare++-tutorial.pdf ..\dynare++\doc\dynare++-tutorial.pdf

 CreateShortcut "${SMLOC}\Documentation.lnk" "$INSTDIR\doc"
SectionEnd

Section "Uninstall"
!insertmacro DETERMINE_CONTEXT

 # First delete the uninstaller
 Delete $INSTDIR\uninstall.exe
 Delete $INSTDIR\README.txt
 Delete $INSTDIR\license.txt
 Rmdir /r $INSTDIR\matlab
 Rmdir /r $INSTDIR\mex
 Rmdir /r $INSTDIR\dynare++
 Rmdir /r $INSTDIR\doc
 # We don't force deletion of installation directory (with /r), to avoid deleting important files
 Rmdir $INSTDIR

 # Delete start menu entries
 Rmdir /r "${SMLOC}"

 # Delete entry in "Add/Remove programs"
 DeleteRegKey SHELL_CONTEXT "${REGLOC}"
SectionEnd

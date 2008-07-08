Dynare version 4
================

For information about how to use Dynare, you should have a look at the
documentation located in the doc/ subdirectory of your Dynare installation (you
should have a shortcut in your Start Menu to access it directly).

Beginners should start with the Dynare userguide. There is also a complete
reference manual documenting all Dynare functions in manual/index.html (however
note that it is a bit outdated).

You can also get more information on the Dynare Wiki at:

  http://www.cepremap.cnrs.fr/DynareWiki

You can also ask questions on the Dynare forums at:

  http://www.dynare.org

Note that Dynare comes with an automated uninstaller, which you can run from
the "Add/Remove Programs" menu of the Control Panel.


Using Dynare with Matlab
------------------------

Dynare requires Matlab version 6.5 or above. With older versions of Matlab, it
may fail or give unexpected results.

To use Dynare, you just have to add the matlab/ subdirectory of your Dynare
installation to Matlab's path. You have two options for doing that:

* Use the addpath command, by typing something like:

    addpath DYNARE_INSTALLATION_DIRECTORY/matlab

  Matlab will not remember this setting next time you run it, and you will have
  to do it again.

* Select the "Set Path" entry in the "File" menu, then click on "Add
  Folder...", and select the matlab/ subdirectory of your Dynare
  installation. Note that you SHOULD NOT use "Add with Subfolders...". Apply
  the settings by clicking on "Save". Note that Matlab will remember this
  setting next time you run it.

You can test your installation by typing "dynare" at the Matlab prompt. This
should give you an error message complaining that you did not specify a MOD
file.


Using Dynare with Octave
------------------------

Dynare is now available for Octave, a free clone of Matlab (see
<http://www.octave.org>).

For installing Octave on your Windows system, go to:

  http://sourceforge.net/project/showfiles.php?group_id=2888&package_id=40078

Then pick the "octave-3.0.1-vs2008-setup.exe" installer. Choose GNUplot
graphical backend (instead of JHandles) during installation. Also note that
this distribution contains a nice text editor, which you can invoke with "edit"
as you would under Matlab.

WARNING: the Octave binary distribution has a bug which makes Octave crash
everytime one types clear all (and therefore everytime one runs Dynare!). A
simple workaround is to type the following command the first time you run
Octave:

  pkg rebuild -noauto ftp ann database

Every time you run Octave, you should type the two following commands:

  addpath DYNARE_INSTALLATION_DIRECTORY/matlab mark_as_command dynare

NOTE: if you don't want to type these two commands every time you run Octave,
you can put them in a file called ".octaverc" in your home directory (generally
"c:\Documents and Settings\USERNAME\"). This file is run by Octave at every
startup.

You can test your installation by typing "dynare" at the Octave prompt. This
should give you an error message complaining that you did not specify a MOD
file.


Copyright notice for software
-----------------------------

Most Dynare source file are copyright "Dynare Team". There are some exceptions
to this, which are described in source file headers when relevant.

Dynare is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program (see gpl.txt). If not, see <http://www.gnu.org/licenses/>.


Copyright notice for documentation
----------------------------------

Permission is granted to copy, distribute and/or modify Dynare documentation
files under the terms of the GNU Free Documentation License, Version 1.2 or any
later version published by the Free Software Foundation; with no Invariant
Sections, no Front-Cover Texts, and no Back-Cover Texts.

A copy of the license can be found in the fdl.txt file.

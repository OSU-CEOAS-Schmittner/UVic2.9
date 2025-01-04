#! /usr/bin/perl -w

#*******************************************************************************
# mk file
#*******************************************************************************

$text1 = "
--------------------------------------------------------------------------------
mk [Command [[Path]Version] [String]]
--------------------------------------------------------------------------------
Command:
 'q' makes executable and submits to a queue (includes mk s)
 'r' makes and runs an executable locally (includes mk e)
 's' makes executable and submit scripts (includes mk e)
 'e' makes an executable (includes mk o) 
 'o' makes all object files (includes mk f)
 'f' makes all compilable files in code
 'file' to make a file (eg. 'mk atmos.o' or 'mk atm.h')
 'c' clears the code directory
 'h' writes help to a log file
--------------------------------------------------------------------------------
Path: model version directory path (default searches \$PATH)
Version: model version directory (only used if not in mk.in or \$HOME/.mk/mk.in)
String: help for mk.in variables containing string (only used if Command = 'h')
--------------------------------------------------------------------------------
Any command will create a local, default mk.in file if one is not found.
--------------------------------------------------------------------------------
";
$text2 = "
General description:
 mk is used to create, compile, link and run code. The program works in a 
 similar way to 'make' in that the minimum amount of code is recompiled 
 when changes are made. mk reads the model configuration from 'mk.in' and 
 'mk.ver' files.  If important definitions are changed, mk will detect these
 and 're-mk' the code.
 
 mk looks for configuration definitions in the initial directory, then the 
 \$HOME/.mk directory and finally in a UVic_ESCM installation subdirectory. Only
 settings in a local mk.in file are used. If a local mk.in file is missing, the
 model version must be specified after the mk command. A default version of 
 mk.in is then copied from the \$HOME/.mk directory and if one is not found 
 there, from a UVic_ESCM installation subdirectory. The Version_Directory 
 variable is always set at the top of a copied default mk.in file so the 
 version only has to be specified once. mk.ver files are treated differently. 
 Settings in all mk.ver files are read but the priority for settings is from the
 mk.ver files in the local directory, then the \$HOME/.mk directory and then a 
 UVic_ESCM installation subdirectory.
  
 mk uses 'Source_Directory' definitions as the basis for all other definitions. 
 Source directories may exist anywhere but usually are under a model
 version directory. Multiple source directories are specified with a bracketed
 index. If a variable does not specify a source directory index, then an index 
 of zero is assumed. Undefined variables will usually get their definitions 
 from lower numbered source directories so specifying a variable without an 
 index specifies that variable for all source directories (until a different 
 definition is specified). The version directory also usually contains a system
 wide 'mk.ver' file which has the system default values. The version directory
 may also contain input files, library files or run files. 
 
 An update level for source directories may be specified and mk will use files 
 in the update level source directories before other specified source 
 directories. The source directories must have the same name under the 
 'Updates_Level' directory for them to be recognized. The updates level source 
 directories must contain all code that has been changed from the original 
 version source directories. No earlier update level sources are used.
 
 mk allows for multiple executables, executable directories, code directories 
 and data directories but only one per source directory. Usually an executable
 directory contains the executable input and output files and code and data
 directories. The code directory contains readable model code and an object and
 dependency directory. It may also have a backup directory. The model is run in
 the executable directory.
 
 You may run executables on the local machine with 'mk r' but you can also
 run multiple executables or queue the executables with 'mk q'. This copies a
 specified perl script run file to 'Executable_File'.run in the executable 
 directory. The perl script run file usually exists under the version directory
 and may be model or system specific. Some mk defined variables are suppled at 
 the beginning of the run file. See examples. The 'Executable_File'.run file 
 should be removed automatically when the executable finishes. If mk finds the
 run file in the executable directory, it will not run another executable until 
 'Executable_File'.run is removed.
 
 The command 'mk h' will print to a log file the commands that can be used, 
 this description and all available preprocessor options.
--------------------------------------------------------------------------------
Operating system derived variables used by mk:
--------------------------------------------------------------------------------
Home:
 A directory used in searching for mk.in files in 'HOME/.mk'. It is derived 
 from the system variable 'HOME'.
Machine_Name:
 This is added to various variables to make them machine dependent. Machine_Name
 is derived from the system variable 'MK_MACHINE_NAME' or if this is not 
 defined, from the command 'uname -n'. Check that this gives the expected 
 response if you wish to use machine dependent variables. Machine_Name is 
 converted to all lower case for use in machine dependent variables. This
 variable may also be specified in mk.in or mk.ver files.
Mk_Path:
 Path that is used to search for Version_Directories. It is derived from the 
 system variable 'PATH'.
Operating_System:
 This is added to various variables to make them operating system dependent.
 Operating_System is derived from the system variable 'MK_OSTYPE' or if this is
 not defined, from from the command 'uname -s'. Check that this gives the 
 expected response if you wish to use operating system dependent variables. 
 Operating_System is converted to all lower case for use in system dependent 
 variables. This variable may also be specified in mk.in or mk.ver files.
User:
 Name used if an email address is requested by a run script. It is derived from
 the system variable 'USER' or read from a mk.in or mk.ver file.
--------------------------------------------------------------------------------
Equality_Character:
--------------------------------------------------------------------------------
 Character that separates a variable name from it's value or setting in a mk.in
 file. It can be changed to any character sequence multiple times in the file.
 Default is '='.
--------------------------------------------------------------------------------
Alphabetical list of variables that may be set for each Source_Directory:
 (see examples in mk.in or mk.ver files)
--------------------------------------------------------------------------------
Auto_Update:
 Flag to automatically update code. Default is 'false'. If true and source code 
 is newer than then the version in code, mk will ask if you want the code 
 updated. Be careful if using this with No_Warnings. 
Backup_Directory:
 Directory that exists under the code directory and contains backup files. 
 Optional. Files are copied from the Code_Directory when the code is cleared 
 with 'mk c'.If the directory does not exist, it will be created.
Change_Mount:
 A list of old and new file system mount points separated by a space. Optional.
 Only used if queueing on a distributed system where a common mount point is
 required but not given by 'pwd'. Directories starting with the old mount point
 are replaced by the new. For example, if 'Change_Mount = /u01 /u', the 
 '/u01/home/' directory would become '/u/home/'. Directories are only changed 
 if the new directory exists.
Code_Directory:
 Directory which contains readable code files. It also usually contains the 
 Dependency_Directory, Object_Directory and Backup Directory. Default is 
 'code'.If the directory does not exist, it will be created.
Code_Extension:
 Code extension after preprocessing (include files retain an Include_Extension).
 Only the first in a list is used. Default is 'Code_Extension_is_undefined'.
Compiler_'Source_Extension' (see Machine_Name and Operating_System):
 Compiler with options. Compiler may be set differently for each source code 
 extension. The last option flag should be set to create object code. The syntax
 for compiling is: 'Compiler_Source_Extension file.Source_Extension'. Compilers 
 may also be set for parallel or specific operating systems or machines.
 Compiler_Parallel definitions are only used if Number_Processors is greater
 than one. Note: if Compiler_Source_Extension is defined, it will overrule other
 definitions. The order of priority is:
  Compiler_'Source_Extension'
  Compiler_Parallel_'Source_Extension'_'Machine_Name'  
  Compiler_Parallel_'Source_Extension'_'Operating_System'
  Compiler_Parallel_'Source_Extension'
  Compiler_Parallel_'Machine_Name'
  Compiler_Parallel_'Operating_System'
  Compiler_Parallel
  Compiler_'Source_Extension_Machine_Name'
  Compiler_'Source_Extension_Operating_System'
  Compiler_'Machine_Name'
  Compiler_'Operating_System'
  Compiler
  Compiler_Default
Compress_Command:
 Command for compressing or taring the code directory. Optional. Syntax is:
 'Compress_Command Code_Directory.Compress_Extension Code_Directory'.
 Code_Directory is removed after compression or tarring.
Compress_Extension:
 Extension for a compressed or tarred Code_Directory.
Data_Directory:
 Directory which will contain data files. Optional. Data files are copied from
 Data_Source if any are missing. If the directory does not exist, it will be 
 created.
Data_Source:
 Directory which will contain data source files. Optional. Default is 
 Version Directory/Data_Directory. If not found locally, mk will look in the 
 Version Directory.
Definitions_Down:
 Flag to spread definitions downward through source directories. Any undefined 
 variable will get its value from from the next highest source directory that
 has the variable defined. Default is true.
Definitions_Up:
 Flag to spread definitions upward through source directories. Default is 
 false.
Dependency_Directory:
 Directory that exists under the code directory and contains dependency files 
 for the mk program. Default is Code_Directory/D.If the directory does not 
 exist, it will be created.
End_String:
 A string that is passed to a run file which can be checked to indicate the
 successful completion of a run. This variable is not passed to the run file if
 set to false or only spaces.
Executable_Directory:
 The directory where the Executable_File will be run. Default is the directory 
 where mk was run.
Executable_File:
 As it sounds. Optional (but recommended!)
Include_Extension:
 Extensions for include files. May be a list separated by spaces. Default is
 Include_Extension_is_undefined.
Input_File:
 A file which is usually read by the Executable_File. Optional.
Libraries (see Machine_Name and Operating_System):
 Libraries which are used in linking code. Optional. The syntax for linking is:
 Linker Executable_File Object_Directory/*.Object_Extension Libraries
 Libraries may also be set for parallel or specific operating systems or 
 machines.Libraries_Parallel definitions are only used if Number_Processors is
 greater than one. Note: if Libraries is defined, it will overrule other 
 definitions. The order of priority is:
  Libraries
  Libraries_Parallel_'Machine_Name'
  Libraries_Parallel_'Operating_System'
  Libraries_Parallel
  Libraries_'Machine_Name'
  Libraries_'Operating_System'
  Libraries_Default
Libraries_Directory (see Machine_Name and Operating_System):
 Library paths which are used in linking code. Optional. These may be a list
 which includes several paths. Libraries_Directory may also be set for parallel
 or specific operating systems or machines. If a directory is not found, mk will 
 look in Version_Directory/lib. Library_Directory_Parallel definitions are only 
 used if Number_Processors is greater than one. Note: if Libraries_Directory is
 defined, it will overrule other definitions. The order of priority is:
  Libraries_Directory
  Libraries_Directory_Parallel_'Machine_Name'
  Libraries_Directory_Parallel_'Operating_System'
  Libraries_Directory_Parallel
  Libraries_Directory_Directory_'Machine_Name'
  Libraries_Directory_'Operating_System'
  Libraries_Directory_Default
Libraries_Directory_Prefix (see Machine_Name and Operating_System):
 Prefix for Libraries_Directory. Default is '-L'. Libraries_Directory_Prefix 
 may also be set for operating systems or machines. Note: if 
 Libraries_Directory_Prefix is defined, it will overrule other definitions. The 
 order of priority is:
  Libraries_Directory_Prefix
  Libraries_Directory_Prefix_'Machine_Name'
  Libraries_Directory_Prefix_'Operating_System'
  Libraries_Directory_Prefix_Default
Linker (see Machine_Name and Operating_System):
 Linker which is used in linking code. Optional. The syntax for linking is: 
 'Linker Executable_File Object_Directory/*.Object_Extension Libraries'
 Linker may also be set for parallel or specific operating systems or machines.
 Linker_Parallel definitions are only used if Number_Processors is greater 
 than one. Note: if Linker is defined, it will overrule other definitions.
 The order of priority is:
  Linker
  Linker_Parallel_'Machine_Name'
  Linker_Parallel_'Operating_System'
  Linker_Parallel
  Linker_'Machine_Name'
  Linker_'Operating_System'
  Linker_Default
Log_File:
 Log file which contains output from mk. Default is log.
Machine_Name:
 This is added to various variables to make them machine dependent. This 
 setting will overrule system definitions.
Mk_Include_File:
 File which contains mk settings for including in code. Optional.
Mk_Include_Start:
 Characters written before settings in Mk_Include_File (eg. print*,). Optional.
Mk_Include_Width:
 Maximum length of settings written to Mk_Include_File. Optional.
Model_Options:
 Model options for the Preprocessor or Compiler. Optional.
Module_Extension (see Machine_Name and Operating_System):
 Extension for modules files that need to be available for compiling. Optional.
 Module_Extension may also be set for operating systems or machines. Note: if 
 Module_Extension is defined, it will overrule other definitions. The order of 
 priority is:
  Module_Extension
  Module_Extension_'Machine_Name'
  Module_Extension_'Operating_System'
  Module_Extension_Default
Nice:
 Nice value sent to the Run_File so the Run_File can set the job priority.
 Optional.
Notify_User:
 Email address that can be used by the run scripts to notify a user when a job
 has finished. If just the domain name is supplied (must start with the 'at' 
 symbol), the user name (see User) will be added. Lists of users will be passed
 on but all spaces are removed. If set to false or only spaces, Notify_User is
 not sent to the Run_File and automatic notification may not be done.
No_Warnings:
 Flag that turns off all warnings from mk. Default is false. If true, code may 
 be deleted without warning. Be careful.
Number_Processors:
 Number of processors. Default is 1. If greater than 1, it will indicate a 
 parallel run.
Operating_System:
 This is added to various variables to make them operating system dependent.
 This setting will overrule system definitions.
Object_Directory:
 Directory that exists under the code directory and contains object files. 
 Default is Code_Directory/O. If the directory does not exist, it will be 
 created.
Object_Extension:
 Extensions for object files. Only the first in a list is used. Default is 
 Object_Extension_is_undefined.
Option_Prefix (see Machine_Name and Operating_System):
 Prefix for Preprocessor options. Default is '-D'. Option_Prefix may also be set
 for operating systems or machines. Note: if Option_Prefix is defined, it will 
 overrule other definitions. The order of priority is:
  Option_Prefix
  Option_Prefix_'Machine_Name'
  Option_Prefix_'Operating_System'
  Option_Prefix_Default
Output_File:
 Standard output from the executable. Optional.
Preprocessor (see Machine_Name and Operating_System):
 Preprocessor used to apply Model_Options. Optional. Preprocessor may also be
 set for operating systems or machines. If Note: if Preprocessor is defined, it 
 will overrule other definitions. The order of priority is:
  Preprocessor
  Preprocessor_'Machine_Name'
  Preprocessor_'Operating_System'
  Preprocessor_Default
Preprocessor_Code:
 Flag to leave code as preprocessed files. Default is false.
Preprocessor_compile:
 Flag to use the compiler for preprocessing. Default is false. If this option
 is used the compiler definitions are modified. mk will look for the string 
 defined by Preprocessor in Compiler_Source_Extension and replace it with 
 Preprocessor followed by all Model_Options (with Option_Prefix included).
 Preprocessor and Option_Prefix may have to be set differently than when not 
 using this option.
Remote_Directory:
 Directory that triggers running on a remote machine. Optional. If
 Executable_Directory starts with Remote_Directory then Machine_Name is set to
 Remote_Machine.
Remote_Machine (Remote_Directory):
 Machine_Name for running on a remote machine. Set to false to never run 
 remotely.
Run_Copy:
 A space separated list of files or directories that are copied to the 
 Run_Directory. Optional.
Run_Directory:
 A remote directory for running the Executable_File. Only used with a Run_File.
 Run_Directory may be used to run the Executable_File on a local scratch 
 directory. Optional.
Run_File (see Machine_Name and Operating_System):
 Files used to submit Executable_File to a queue. Optional. Run_File may also
 be set for operating systems or machines. Run_File_Parallel definitions are 
 only used if Number_Processors is greater than one. If Note: if Run_File is
 defined, it will overrule other definitions. The order of priority is:
  Run_File
  Run_File_Parallel_'Machine_Name'
  Run_File_Parallel_'Operating_System'
  Run_File_Parallel
  Run_File_'Machine_Name'
  Run_File_'Operating_System'
  Run_File_Default
Run_File_Local (see Run_File):
 Run_File to be used if running locally (with 'mk r'). Set to false to run
 locally without a Run_File.
Search_Path:
 Path that mk will search for a Version_Directory. Optional. May be a space 
 separated list.
Source_Comment_Character:
 Unless set to false, mk will add the Source_Comment_Character followed by the
 name of the source file to the beginning of code and include files. Optional.
Source_Directory:
 Directories which contain source or include files. If not found locally, mk
 will look in the Version Directory.
Source_Extension:
 Extensions for source files. May be a list separated by spaces. Default is
 Source_Extension_is_undefined.
Step_Start_File:
 Unless undefined or set to false, an executable that will be run before each
 step in a multi-step run. Optional.
Step_End_File:
 Unless undefined or set to false, an executable that will be run after each
 step in a multi-step run. Optional.
Uncompress_Command:
 Command for uncompressing or untaring the code directory. Syntax is:
 Compress_Command Code_Directory.Compress_Extension
Updates_Level:
 Directory with updated source code. Optional. Default is latest. Often set
 to a number which determines the updates level. If the directory is not found
 under 'updates' locally it will be searched for in Version_Directory/updates. 
 If set to latest, Version_Directory/updates will be searched for the highest 
 number. 
User
 Name used if an email address is requested by a run script. This setting will
 overrule system definitions.
Version_Directory:
 Directory which determines the model version. Must be defined but may be set 
 as a command option.
--------------------------------------------------------------------------------
";
$line = "
--------------------------------------------------------------------------------
";
$text1 =~ s/^\n//; $text2 =~ s/^\n//; $line =~ s/^\n//;

#-------------------------------------------------------------------------------
# set some defaults
#-------------------------------------------------------------------------------

$Code_Directory = "code";
$Object_Directory = "O";
$Dependency_Directory = "D";
$Log_File = "mk.log";
$Executable_File = "mk.exe";
$Updates_Level = "latest";
$Definitions_Up = "true";
$Definitions_Down = "false";
$Auto_Update = "false";
$No_Warnings = "false";
$Preprocessor_Code = "false";
$Preprocessor_Compile = "false";
$Notify_User = "false";
$End_String = "false";
$Equality_Character = "=";
$Libraries_Directory_Prefix = "-L";
$Option_Prefix = "-D";

#-------------------------------------------------------------------------------
# get relevant environment variables
#-------------------------------------------------------------------------------

$Home = $ENV{"HOME"};
unless ( $Home ) { $Home = " " }
$Operating_System = $ENV{"MK_OPERATING_SYSTEM"};
unless ( $Operating_System ) { $Operating_System = `uname -s` }
unless ( $Operating_System ) { $Operating_System = " " }
chomp $Operating_System;
$Operating_System = lc $Operating_System;
$Machine_Name = $ENV{"MK_MACHINE_NAME"};
unless ( $Machine_Name ) { $Machine_Name = `uname -n` }
unless ( $Machine_Name ) { $Machine_Name = " " }
chomp $Machine_Name;
$Machine_Name = lc $Machine_Name;
$Mk_Path = $ENV{"PATH"};
unless ( $Mk_Path ) { $Mk_Path = " " }
$User = $ENV{"USER"};
unless ( $User ) { $User = " " }

use File::Copy;
use Cwd;
$Initial_Directory = cwd();
$Cmd = $ARGV[0];
$Version = $ARGV[1];
$String = $ARGV[2];
unless ( $String ) { $String = $ARGV[1] }

#-------------------------------------------------------------------------------
# read mk.in input file (look in local, home and version directories)
#-------------------------------------------------------------------------------

if ( $Version ) {
 if ( ! -e  $Version ) {
  $Version = find_path ( $Version, " ", $Mk_Path );
 }
}

# if local mk.in is missing get one from $HOME/.mk
if ( ! -e "mk.in" ) {
 if ( -e "$Home/.mk/mk.in"  && $Version ) {
  undef @out; 
  push ( @out, "Version_Directory $Equality_Character $Version\n" );
  open IN, "$Home/.mk/mk.in";
  while (<IN>) { push ( @out, "$_" ) };
  close IN;
  open OUT, ">mk.in"; foreach (@out) { print OUT "$_" }; close OUT;
  read_mk ( "mk.in" ); %Mk = %MkH;
 }
}

# if local mk.in is still missing get one from $Version or $Version/run
if ( ! -e "mk.in" && $Version ) {
 if ( -e  "$Version/mk.in" ) {
  undef @out; 
  push ( @out, "Version_Directory $Equality_Character $Version\n" );
  open IN, "$Version/mk.in";
  while (<IN>) { push ( @out, "$_" ) };
  close IN;
  open OUT, ">mk.in"; foreach (@out) { print OUT "$_" }; close OUT;
  read_mk ( "mk.in" ); %Mk = %MkH;
 }
 else {
  if ( -e "$Version/run/mk.in" ) {
   undef @out; 
   push ( @out, "Version_Directory $Equality_Character $Version\n" );
   open IN, "$Version/run/mk.in";
   while (<IN>) { push ( @out, "$_" ) };
   close IN;
   open OUT, ">mk.in"; foreach (@out) { print OUT "$_" }; close OUT;
   read_mk ( "mk.in" ); %Mk = %MkH;
  }
 }
}
read_mk ( "mk.in" ); %Mk = %MkH;

#-------------------------------------------------------------------------------
# read mk.ver input file (look in local, home and version directories)
#-------------------------------------------------------------------------------

read_mk ( "mk.ver" ); %Mk = %MkH;
read_mk ( "$Home/.mk/mk.ver" ); %Mk = %MkH;

# search for input file in $Version if it is set as an argument
if ( $Version ) {
 if ( ! -e  $Version ) {
# search through all possible search paths
  for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
   $Version = find_path ( $Version, $Mk{Search_Path}[$idx], $Mk_Path );
   if ( -e  "$Version/mk.ver" ) {
    read_mk ( "$Version/mk.ver", %Mk ); %Mk = %MkH;
   }
   else {
    if ( -e  "$Version/run/mk.ver" ) {
     read_mk ( "$Version/run/mk.ver", %Mk ); %Mk = %MkH;
    }
   }
  }
 }
 else {
  if( -e  "$Version/mk.ver" ) {
   read_mk ( "$Version/mk.ver", %Mk ); %Mk = %MkH;
  }
  else {
   if( -e  "$Version/run/mk.ver" ) {
    read_mk ( "$Version/run/mk.ver", %Mk ); %Mk = %MkH;
   }
  }
 }
}

# search for input files in any Version_Directory
for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
 if ( $Mk{Version_Directory}[$idx] ) {
  $Version = $Mk{Version_Directory}[$idx];
  if ( ! -e  $Version) {
   $Version = find_path ( $Version, $Mk{Search_Path}[$idx], $Mk_Path );
  }
  if ( -e  "$Version/mk.ver" ) {
   read_mk ( "$Version/mk.ver", %Mk ); %Mk = %MkH;
  }
  else {
   if ( -e  "$Version/run/mk.ver" ) {
    read_mk ( "$Version/run/mk.ver", %Mk ); %Mk = %MkH;
   }
  }
 }
}

# skip if no command argument
if ( $Cmd ) {
 $error = "false";
 if ( $Version ) { if ( ! -e $Version ) { $error = "true" } }
 else { $error = "true" }
# stop if $Version is not found
 if ( $error eq "true" ) {
  print "$text1\n";
  print "Mk can not find a Version_Directory.\n\n";
  print "  You may need to add the model version after the mk command,\n";
  print "  set your PATH or add a Version_Directory definition to mk.in.\n";
  exit;
 }
}

#-------------------------------------------------------------------------------
# Machine_Name, Operating_System and User can be set in mk.in or mk.ver but 
# only the first definition is used.
#-------------------------------------------------------------------------------

 if ( $Mk{Machine_Name}[0] ) { $Machine_Name = $Mk{Machine_Name}[0] }
 $Machine_Name = lc $Machine_Name;
 if ( $Mk{Operating_System}[0] ) { $Operating_System = $Mk{Operating_System}[0] }
 $Operating_System = lc $Operating_System;
 if ( $Mk{User}[0] ) { $User = $Mk{User}[0] }

#-------------------------------------------------------------------------------
# determine settings for all source directories 
#-------------------------------------------------------------------------------

for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {


#-------------------------------------------------------------------------------
# search for and set the Version_Directory
#-------------------------------------------------------------------------------

 if ( $Mk{Version_Directory}[$idx] ) {
  $Version = $Mk{Version_Directory}[$idx];
  if ( ! -e  $Version) {
   $Version = find_path ( $Version, $Mk{Search_Path}[$idx], $Mk_Path );
  }
  if ( -e  $Version) { $Mk{Version_Directory}[$idx] = $Version }
  unless ( -e $Mk{Version_Directory}[$idx] ) {
   print "Mk can not find a Version_Directory for Source_Directory $idx.\n";
   exit;
  }
 }
 elsif ( $Version ) { $Mk{Version_Directory}[$idx] = $Version }
 unless ( $Mk{Version_Directory}[$idx] ) {
  print "Mk can not find a Version_Directory for Source_Directory $idx.\n"; exit
 }

#-------------------------------------------------------------------------------
# set various files and directories
#-------------------------------------------------------------------------------

# executable directory
 unless ( $Mk{Executable_Directory}[$idx] ) { 
  $Mk{Executable_Directory}[$idx] = $Initial_Directory;
 }
 unless ( -e $Mk{Executable_Directory}[$idx] ) { 
  $Mk{Executable_Directory}[$idx] = $Initial_Directory;
 }
 $Mk{Executable_Directory}[$idx] = set_var ( $Mk{Executable_Directory}[$idx],
                                             $Initial_Directory );
# executable file
 unless ( $Mk{Executable_File}[$idx] ) { 
  $Mk{Executable_File}[$idx] = $Executable_File;
 }
 $Mk{Executable_File}[$idx] = set_var ( $Mk{Executable_File}[$idx],
                                        $Mk{Executable_Directory}[$idx] );

# check for remote machine run
 if ( $Mk{Remote_Directory}[$idx] && $Mk{Remote_Machine}[$idx] ) {
  $a = index ($Mk{Executable_Directory}[$idx], $Mk{Remote_Directory}[$idx] );
  if ( $a == 0 && $Mk{Remote_Machine}[$idx] ne "false") {
   $Machine_Name = lc $Mk{Remote_Machine}[$idx];
  }
 }

# code directory
 unless ( $Mk{Code_Directory}[$idx] ) {
  $Mk{Code_Directory}[$idx] = $Code_Directory;
 }
 $Mk{Code_Directory}[$idx] = set_var ( $Mk{Code_Directory}[$idx], 
                                       $Mk{Executable_Directory}[$idx] );
# log file
 unless ( $Mk{Log_File}[$idx] ) { $Mk{Log_File}[$idx] = $Log_File }
 $Mk{Log_File}[$idx] = set_var ( $Mk{Log_File}[$idx],
                                 $Mk{Executable_Directory}[$idx] );
# input file
 $Mk{Input_File}[$idx] = set_var ( $Mk{Input_File}[$idx], 
                                   $Mk{Executable_Directory}[$idx] );
# output file
 $a = set_var ( $Mk{Output_File}[$idx], $Mk{Executable_Directory}[$idx] );
 $Mk{Output_File}[$idx] = set_var ( $Mk{Output_File}[$idx],
                                    $Mk{Executable_Directory}[$idx] );
# data source directory
 if ( $Mk{Data_Directory}[$idx] && ! $Mk{Data_Source}[$idx] ) {
  $a = set_var ( $Mk{Data_Directory}[$idx], $Mk{Version_Directory}[$idx] );
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $Mk{Data_Directory}[$idx], $Mk{Version_Directory}[$idx] );
  }
  if ( -e $a ) { $Mk{Data_Source}[$idx] = $a }
 }
 if ( $Mk{Data_Source}[$idx] ) {
  if ( ! -e $Mk{Data_Source}[$idx] ) {
   $Mk{Data_Source}[$idx] = set_var ( $Mk{Data_Source}[$idx],
                                      $Mk{Version_Directory}[$idx] );
  }
 }
 
# data directory
 $Mk{Data_Directory}[$idx] = set_var ( $Mk{Data_Directory}[$idx],
                                       $Mk{Executable_Directory}[$idx] );
# email
  unless ( $Mk{Notify_User}[$idx] ) { 
   $Mk{Notify_User}[$idx] = $Notify_User
  }
  if ( $Mk{Notify_User}[$idx] ) {
   if ( index ( $Mk{Notify_User}[$idx], "@" ) == 0 ) {
     $Mk{Notify_User}[$idx] = "$User$Mk{Notify_User}[$idx]";
   }
   $Mk{Notify_User}[$idx] =~ s/ //g
  }
  if ( $Mk{Notify_User}[$idx] eq "false" ) { undef $Mk{Notify_User}[$idx] }
  
# end string
  unless ( $Mk{End_String}[$idx] ) { $Mk{End_String}[$idx] = $End_String }
  if ( $Mk{End_String}[$idx] eq "false" ) { undef $Mk{End_String}[$idx] }

# nice value
  if ( $Mk{Nice}[$idx] ) {
   if ( $Mk{Nice}[$idx] eq "false" ) { undef $Mk{Nice}[$idx] }
  }
  
# number of processors
 if ( $Mk{Number_Processors}[$idx] ) {
  if ( $Mk{Number_Processors}[$idx] <= 1 ) {
   undef $Mk{Number_Processors}[$idx];
  }
 }

#-------------------------------------------------------------------------------
# look for a "run" file
#-------------------------------------------------------------------------------

 if ( $Cmd ) {
  if ($Cmd eq "e" || $Cmd eq "r" || $Cmd eq "s" || $Cmd eq "q" ) {
   if ( $Mk{Executable_File}[$idx] ) {
    if ( -e "$Mk{Executable_File}[$idx].run" ) {
     print "Mk found a run file. A model may be running in this directory.\n";
     print "If a model is not running, then resubmit after deleting\n";
     print "$Mk{Executable_File}[$idx].run\n";
     exit;
    }
   }
  }
 }

#-------------------------------------------------------------------------------
# look for Step_Start_File and Step_End_File
#-------------------------------------------------------------------------------

 if ( $Mk{Step_Start_File}[$idx] ) {
  $a = $Mk{Step_Start_File}[$idx];
  if ( ! -e $a ) {
   $a = set_var ( $Mk{Step_Start_File}[$idx], $Mk{Executable_Directory}[$idx] );
  }
  if ( ! -e $a ) {
   $a = set_var ( $Mk{Step_Start_File}[$idx], "$Home/.mk" );
  }
  if ( ! -e $a ) {
   $a = set_var ( $Mk{Step_Start_File}[$idx],
                  "$Mk{Version_Directory}[$idx]/run" );
  }
  if ( ! -e $a ) {
   $a = set_var ( $Mk{Step_Start_File}[$idx], $Mk{Version_Directory}[$idx] );
  }
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $a, $Mk{Executable_Directory}[$idx] );
  }
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $a, "$Home/.mk" );
  }
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $a, "$Mk{Version_Directory}[$idx]/run" );
  }
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $a, $Mk{Version_Directory}[$idx] );
  }
  if ( $Mk{Step_Start_File}[$idx] eq "false" ) { 
   undef $Mk{Step_Start_File}[$idx] 
  }
  else {$Mk{Step_Start_File}[$idx] = $a }
 }
 
 if ( $Mk{Step_End_File}[$idx] ) {
  $a = $Mk{Step_End_File}[$idx];
  if ( ! -e $a ) {
   $a = set_var ( $Mk{Step_End_File}[$idx], $Mk{Executable_Directory}[$idx] );
  }
  if ( ! -e $a ) {
   $a = set_var ( $Mk{Step_End_File}[$idx], "$Home/.mk" );
  }
  if ( ! -e $a ) {
   $a = set_var ( $Mk{Step_End_File}[$idx],
                  "$Mk{Version_Directory}[$idx]/run" );
  }
  if ( ! -e $a ) {
   $a = set_var ( $Mk{Step_End_File}[$idx], $Mk{Version_Directory}[$idx] );
  }
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $a, $Mk{Executable_Directory}[$idx] );
  }
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $a, "$Home/.mk" );
  }
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $a, "$Mk{Version_Directory}[$idx]/run" );
  }
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $a, $Mk{Version_Directory}[$idx] );
  }
  if ( $Mk{Step_End_File}[$idx] eq "false" ) {
   undef $Mk{Step_End_File}[$idx]
  }
  else {$Mk{Step_End_File}[$idx] = $a }
 }

#-------------------------------------------------------------------------------
# set directories which can only exist below the code directory
#-------------------------------------------------------------------------------

 if ( $Mk{Backup_Directory}[$idx] ) {
  $a = $Mk{Backup_Directory}[$idx]; $a =~ s/.*\///g; $a =~ s/ //g;
  $Mk{Backup_Directory}[$idx] = "$Mk{Code_Directory}[$idx]/$a"
 }
 
 unless ($Mk{Dependency_Directory}[$idx] ) {
  $Mk{Dependency_Directory}[$idx] = $Dependency_Directory;
 }
 $a = $Mk{Dependency_Directory}[$idx]; $a =~ s/.*\///g; $a =~ s/ //g;
 $Mk{Dependency_Directory}[$idx] = "$Mk{Code_Directory}[$idx]/$a";

 unless ($Mk{Object_Directory}[$idx] ) {
  $Mk{Object_Directory}[$idx] = $Object_Directory;
 }
 $a = $Mk{Object_Directory}[$idx]; $a =~ s/.*\///g; $a =~ s/ //g;
 $Mk{Object_Directory}[$idx] = "$Mk{Code_Directory}[$idx]/$a";
 
#-------------------------------------------------------------------------------
# if a source is not found, try initial, executable and version directories
#-------------------------------------------------------------------------------

 $a = set_var ( $Mk{Source_Directory}[$idx], $Initial_Directory );
 if ( ! -e $a ) {
  $a = set_var ( $Mk{Source_Directory}[$idx], $Mk{Executable_Directory}[$idx] );
 }
 if ( ! -e $a ) {
  $a = set_var ( $Mk{Source_Directory}[$idx], $Mk{Version_Directory}[$idx] );
 }
 $Mk{Source_Directory}[$idx] = $a;
 
#-------------------------------------------------------------------------------
# find current Updates_Level
#-------------------------------------------------------------------------------

 unless ( $Mk{Updates_Level}[$idx] ) { 
  $Mk{Updates_Level}[$idx] = $Updates_Level
 }
 if ( $Mk{Updates_Level}[$idx] eq $Updates_Level ) {
  @files = sort ( <$Mk{Version_Directory}[$idx]/updates/*> );
  if (@files) { $Mk{Updates_Level}[$idx] = @files[ $#{ @files } ] }
 }
 $a = "$Initial_Directory/updates/";
 $a = set_var ( $Mk{Updates_Level}[$idx], $a );
 if ( ! -e $a ) {
  $a = "Executable_Directory}[$idx]/updates/";
  $a = set_var ( $Mk{Updates_Level}[$idx], $a );
 }
 if ( ! -e $a ) {
  $a = "$Mk{Version_Directory}[$idx]/updates/";
  $a = set_var ( $Mk{Updates_Level}[$idx], $a );
 }
 $Mk{Updates_Level}[$idx] = $a;
 if ( $Mk{Updates_Level}[$idx] eq $Updates_Level ) {
  $Mk{Updates_Level}[$idx] = " ";
 }
 
#-------------------------------------------------------------------------------
# sort model options
#-------------------------------------------------------------------------------

 unless ( $Mk{Model_Options}[$idx] ) { $Mk{Model_Options}[$idx] = " " }
 @options = sort ( split ( / /, $Mk{Model_Options}[$idx] ) );
 $Mk{Model_Options}[$idx] = join ( ' ', @options );

#-------------------------------------------------------------------------------
# set option_prefix
#-------------------------------------------------------------------------------

 $key = "Option_Prefix";
 $key_mn = "Option_Prefix_$Machine_Name";
 $key_os = "Option_Prefix_$Operating_System";
 $key_d = "Option_Prefix_Default";
 undef $a;
 foreach $k ( keys %Mk ) { 
  if ( $k eq $key) { unless ( $a ) { $a = $Mk{$key}[$idx] } } 
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_mn ) { unless ( $a ) { $a = $Mk{$key_mn}[$idx] } } 
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_os) { unless ( $a ) { $a = $Mk{$key_os}[$idx] } }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_d) { unless ( $a ) { $a = $Mk{$key_d}[$idx] } }
 }
 if ( $a ) { $Mk{$key}[$idx] = $a }

 unless ( $Mk{Option_Prefix}[$idx] ) { 
  $Mk{Option_Prefix}[$idx] = $Option_Prefix;
 }
 
#-------------------------------------------------------------------------------
# set preprocessor
#-------------------------------------------------------------------------------

 $key = "Preprocessor";
 $key_mn = "Preprocessor_$Machine_Name";
 $key_os = "Preprocessor_$Operating_System";
 $key_d = "Preprocessor_Default";
 undef $a;
 foreach $k ( keys %Mk ) { 
  if ( $k eq $key) { unless ( $a ) { $a = $Mk{$key}[$idx] } } 
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_mn ) { unless ( $a ) { $a = $Mk{$key_mn}[$idx] } } 
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_os) { unless ( $a ) { $a = $Mk{$key_os}[$idx] } }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_d) { unless ( $a ) { $a = $Mk{$key_d}[$idx] } }
 }
 if ( $a ) { $Mk{$key}[$idx] = $a }
 
#-------------------------------------------------------------------------------
# set module_extension
#-------------------------------------------------------------------------------

 $key = "Module_Extension";
 $key_mn = "Module_Extension_$Machine_Name";
 $key_os = "Module_Extension_$Operating_System";
 $key_d = "Module_Extension_Default";
 undef $a;
 foreach $k ( keys %Mk ) { 
  if ( $k eq $key) { unless ( $a ) { $a = $Mk{$key}[$idx] } } 
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_mn ) { unless ( $a ) { $a = $Mk{$key_mn}[$idx] } } 
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_os) { unless ( $a ) { $a = $Mk{$key_os}[$idx] } }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_d) { unless ( $a ) { $a = $Mk{$key_d}[$idx] } }
 }
 if ( $a ) { $Mk{$key}[$idx] = $a }
 
#-------------------------------------------------------------------------------
# if not defined, define extensions as undefined.
#-------------------------------------------------------------------------------

 unless ( $Mk{Source_Extension}[$idx] ) {
  $Mk{Source_Extension}[$idx] = "Source_Extension_is_undefined";
 }
 unless ( $Mk{Include_Extension}[$idx] ) {
  $Mk{Include_Extension}[$idx] = "Include_Extension_is_undefined";
 }
 unless ( $Mk{Code_Extension}[$idx] ) {
  $Mk{Code_Extension}[$idx] = "Code_Extension_is_undefined";
 }
 unless ( $Mk{Object_Extension}[$idx] ) {
  $Mk{Object_Extension}[$idx] = "Object_Extension_is_undefined";
 }
 unless ( $Mk{Module_Extension}[$idx] ) {
  $Mk{Module_Extension}[$idx] = "Module_Extension_is_undefined"; 
 }
 
#-------------------------------------------------------------------------------
# set run_file
#-------------------------------------------------------------------------------

 $key = "Run_File";
 $key_p_mn = "Run_File_Parallel_$Machine_Name";
 $key_p_os = "Run_File_Parallel_$Operating_System";
 $key_p = "Run_File_Parallel";
 $key_mn = "Run_File_$Machine_Name";
 $key_os = "Run_File_$Operating_System";
 $key_d = "Run_File_Default";
 undef $a;
 foreach $k ( keys %Mk ) {
  if ( $k eq $key) { unless ( $a ) { $a = $Mk{$key}[$idx] } }
 }
 if ( $Mk{Number_Processors}[$idx] ) {
  if ( $Mk{Number_Processors}[$idx] > 1 ) { 
   foreach $k ( keys %Mk ) {
    if ( $k eq $key_p_mn) { unless ( $a ) { $a = $Mk{$key_p_mn}[$idx] } }
   }
   foreach $k ( keys %Mk ) {
    if ( $k eq $key_p_os) { unless ( $a ) { $a = $Mk{$key_p_os}[$idx] } }
   }
   foreach $k ( keys %Mk ) { 
    if ( $k eq $key_p) { unless ( $a ) { $a = $Mk{$key_p}[$idx] } } 
   }
  }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_mn) { unless ( $a ) { $a = $Mk{$key_mn}[$idx] } }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_os) { unless ( $a ) { $a = $Mk{$key_os}[$idx] } }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_d) { unless ( $a ) { $a = $Mk{$key_d}[$idx] } }
 }
 if ( $a ) { $Mk{$key}[$idx] = $a }

 if ( $Mk{Run_File}[$idx] ) {
  $a = $Mk{Run_File}[$idx];
  if ( ! -e $a ) {
   $a = set_var ( $Mk{Run_File}[$idx], $Mk{Executable_Directory}[$idx] );
  }
  if ( ! -e $a ) {
   $a = set_var ( $Mk{Run_File}[$idx], "$Home/.mk" );
  }
  if ( ! -e $a ) {
   $a = set_var ( $Mk{Run_File}[$idx], "$Mk{Version_Directory}[$idx]/run" );
  }
  if ( ! -e $a ) {
   $a = set_var ( $a, $Mk{Version_Directory}[$idx] );
  }
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $a, $Mk{Executable_Directory}[$idx] );
  }
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $a, "$Home/.mk" );
  }
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $a, "$Mk{Version_Directory}[$idx]/run" );
  }
  if ( ! -e $a ) {
   $a =~ s/.*\///g;
   $a = set_var ( $a, $Mk{Version_Directory}[$idx] );
  }
  if ( -e $a ) { $Mk{Run_File}[$idx] = $a }
 }
 
#-------------------------------------------------------------------------------
# set libraries
#-------------------------------------------------------------------------------

 $key = "Libraries";
 $key_p_mn = "Libraries_Parallel_$Machine_Name";
 $key_p_os = "Libraries_Parallel_$Operating_System";
 $key_p = "Libraries_Parallel";
 $key_mn = "Libraries_$Machine_Name";
 $key_os = "Libraries_$Operating_System";
 $key_d = "Libraries_Default";
 undef $a;
 foreach $k ( keys %Mk ) { 
  if ( $k eq $key) { unless ( $a ) { $a = $Mk{$key}[$idx] } } 
 }
 if ( $Mk{Number_Processors}[$idx] ) {
  if ( $Mk{Number_Processors}[$idx] > 1 ) { 
   foreach $k ( keys %Mk ) {
    if ( $k eq $key_p_mn) { unless ( $a ) { $a = $Mk{$key_p_mn}[$idx] } }
   }
   foreach $k ( keys %Mk ) {
    if ( $k eq $key_p_os) { unless ( $a ) { $a = $Mk{$key_p_os}[$idx] } }
   }
   foreach $k ( keys %Mk ) { 
    if ( $k eq $key_p) { unless ( $a ) { $a = $Mk{$key_p}[$idx] } } 
   }
  }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_mn) { unless ( $a ) { $a = $Mk{$key_mn}[$idx] } }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_os) { unless ( $a ) { $a = $Mk{$key_os}[$idx] } }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_d) { unless ( $a ) { $a = $Mk{$key_d}[$idx] } }
 }
 if ( $a ) { $Mk{$key}[$idx] = $a }
 
#-------------------------------------------------------------------------------
# set libraries directory
#-------------------------------------------------------------------------------

 $key = "Libraries_Directory";
 $key_p_mn = "Libraries_Directory_Parallel_$Machine_Name";
 $key_p_os = "Libraries_Directory_Parallel_$Operating_System";
 $key_p = "Libraries_Directory_Parallel";
 $key_mn = "Libraries_Directory_$Machine_Name";
 $key_os = "Libraries_Directory_$Operating_System";
 $key_d = "Libraries_Directory_Default";
 undef $a;
 foreach $k ( keys %Mk ) { 
  if ( $k eq $key) { unless ( $a ) { $a = $Mk{$key}[$idx] } } 
 }
 if ( $Mk{Number_Processors}[$idx] ) {
  if ( $Mk{Number_Processors}[$idx] > 1 ) { 
   foreach $k ( keys %Mk ) {
    if ( $k eq $key_p_mn) { unless ( $a ) { $a = $Mk{$key_p_mn}[$idx] } }
   }
   foreach $k ( keys %Mk ) {
    if ( $k eq $key_p_os) { unless ( $a ) { $a = $Mk{$key_p_os}[$idx] } }
   }
   foreach $k ( keys %Mk ) { 
    if ( $k eq $key_p) { unless ( $a ) { $a = $Mk{$key_p}[$idx] } } 
   }
  }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_mn) { unless ( $a ) { $a = $Mk{$key_mn}[$idx] } }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_os) { unless ( $a ) { $a = $Mk{$key_os}[$idx] } }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_d) { unless ( $a ) { $a = $Mk{$key_d}[$idx] } }
 }

 if ( $a ) { $Mk{$key}[$idx] = $a }
 else { $Mk{$key}[$idx] = "$Mk{Version_Directory}[$idx]/lib" }

 @list = split ( / /,$Mk{$key}[$idx] );
 undef $a;
 foreach (@list) {
  if ( ! -e $_ ) {
   $list = $_;
   $_ = set_var ( $list,"$Mk{Version_Directory}[$idx]" );
   if ( ! -e $_ ) { $_ = set_var ( $list,"$Mk{Version_Directory}[$idx]/lib" ) }
   if ( -e $_ ) { if ( $a ) { $a = "$a $_" } else { $a = "$_" } }
  }
  else { if ( $a ) { $a = "$a $_" } else { $a = "$_" } }
 }
 $Mk{$key}[$idx] = $a;

#-------------------------------------------------------------------------------
# set libraries directory prefix
#-------------------------------------------------------------------------------

 $key = "Libraries_Directory_Prefix";
 $key_mn = "Libraries_Directory_Prefix_$Machine_Name";
 $key_os = "Libraries_Directory_Prefix_$Operating_System";
 $key_d = "Libraries_Directory_Prefix_Default";
 undef $a;
 foreach $k ( keys %Mk ) { 
  if ( $k eq $key) { unless ( $a ) { $a = $Mk{$key}[$idx] } } 
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_mn ) { unless ( $a ) { $a = $Mk{$key_mn}[$idx] } } 
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_os) { unless ( $a ) { $a = $Mk{$key_os}[$idx] } }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_d) { unless ( $a ) { $a = $Mk{$key_d}[$idx] } }
 }
 if ( $a ) { $Mk{$key}[$idx] = $a }

 unless ( $Mk{Libraries_Directory_Prefix}[$idx] ) { 
  $Mk{Libraries_Directory_Prefix}[$idx] = $Libraries_Directory_Prefix;
 }
#-------------------------------------------------------------------------------
# set linker
#-------------------------------------------------------------------------------

 $key = "Linker";
 $key_p_mn = "Linker_Parallel_$Machine_Name";
 $key_p_os = "Linker_Parallel_$Operating_System";
 $key_p = "Linker_Parallel";
 $key_mn = "Linker_$Machine_Name";
 $key_os = "Linker_$Operating_System";
 $key_d = "Linker_Default";
 undef $a;
 foreach $k ( keys %Mk ) {
  if ( $k eq $key) { unless ( $a ) { $a = $Mk{$key}[$idx] } }
 }
 if ( $Mk{Number_Processors}[$idx] ) {
  if ( $Mk{Number_Processors}[$idx] > 1 ) { 
   foreach $k ( keys %Mk ) {
    if ( $k eq $key_p_mn) { unless ( $a ) { $a = $Mk{$key_p_mn}[$idx] } } 
   }
   foreach $k ( keys %Mk ) {
    if ( $k eq $key_p_os) { unless ( $a ) { $a = $Mk{$key_p_os}[$idx] } } 
   }
   foreach $k ( keys %Mk ) {
    if ( $k eq $key_p) { unless ( $a ) { $a = $Mk{$key_p}[$idx] } }
   }
  }
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_mn) { unless ( $a ) { $a = $Mk{$key_mn}[$idx] } } 
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_os) { unless ( $a ) { $a = $Mk{$key_os}[$idx] } } 
 }
 foreach $k ( keys %Mk ) {
  if ( $k eq $key_d) { unless ( $a ) { $a = $Mk{$key_d}[$idx] } } 
 }
 if ( $a ) { $Mk{$key}[$idx] = $a }
 
#-------------------------------------------------------------------------------
# set compiler
#-------------------------------------------------------------------------------

 @extension_list = split(/ /,$Mk{Source_Extension}[$idx]);
# set compiler for each Source_Extension
 foreach $extension (@extension_list) {
  $extension =~ s/ +//g;
  $key_ex = "Compiler_$extension";
  $key_p_ex_mn = "Compiler_Parallel_$extension"."_$Machine_Name";
  $key_p_ex_os = "Compiler_Parallel_$extension"."_$Operating_System";  
  $key_p_ex = "Compiler_Parallel_$extension";
  $key_p_mn = "Compiler_Parallel_$Machine_Name";
  $key_p_os = "Compiler_Parallel_$Operating_System";
  $key_p = "Compiler_Parallel";
  $key_ex_mn = "Compiler_$extension"."_$Machine_Name";
  $key_ex_os = "Compiler_$extension"."_$Operating_System";
  $key_mn = "Compiler_$Machine_Name";
  $key_os = "Compiler_$Operating_System";
  $key = "Compiler";
  $key_d_ex = "Compiler_Default_$extension";
  undef $a;
  foreach $k ( keys %Mk ) { 
   if ( $k eq $key_ex) { unless ( $a ) { $a = $Mk{$key_ex}[$idx] } }
  }
  if ( $Mk{Number_Processors}[$idx] ) {
   if ( $Mk{Number_Processors}[$idx] > 1 ) { 
    foreach $k ( keys %Mk ) { 
     if ( $k eq $key_ex) { unless ( $a ) { $a = $Mk{$key_ex}[$idx] } }
    }
    foreach $k ( keys %Mk ) {
     if ( $k eq $key_p_ex_mn) { unless ( $a ) { $a = $Mk{$key_p_ex_mn}[$idx] } }
    }
    foreach $k ( keys %Mk ) {
     if ( $k eq $key_p_ex_os) { unless ( $a ) { $a = $Mk{$key_p_ex_os}[$idx] } }
    }
    foreach $k ( keys %Mk ) { 
     if ( $k eq $key_p_ex) { unless ( $a ) { $a = $Mk{$key_p_ex}[$idx] } }
    }
    foreach $k ( keys %Mk ) { 
     if ( $k eq $key_p_mn) { unless ( $a ) { $a = $Mk{$key_p_mn}[$idx] } }
    }
    foreach $k ( keys %Mk ) { 
     if ( $k eq $key_p_os) { unless ( $a ) { $a = $Mk{$key_p_os}[$idx] } }
    }
    foreach $k ( keys %Mk ) { 
     if ( $k eq $key_p) { unless ( $a ) { $a = $Mk{$key_p}[$idx] } }
    }
   }
  }
  foreach $k ( keys %Mk ) { 
   if ( $k eq $key_ex) { unless ( $a ) { $a = $Mk{$key_ex}[$idx] } }
  }
  foreach $k ( keys %Mk ) {
   if ( $k eq $key_ex_mn) { unless ( $a ) { $a = $Mk{$key_ex_mn}[$idx] } }
  }
  foreach $k ( keys %Mk ) {
   if ( $k eq $key_ex_os) { unless ( $a ) { $a = $Mk{$key_ex_os}[$idx] } }
  }
  foreach $k ( keys %Mk ) { 
   if ( $k eq $key_mn) { unless ( $a ) { $a = $Mk{$key_mn}[$idx] } }
  }
  foreach $k ( keys %Mk ) { 
   if ( $k eq $key_os) { unless ( $a ) { $a = $Mk{$key_os}[$idx] } }
  }
  foreach $k ( keys %Mk ) { 
   if ( $k eq $key) { unless ( $a ) { $a = $Mk{$key}[$idx] } }
  }
  foreach $k ( keys %Mk ) { 
   if ( $k eq $key_d_ex) { unless ( $a ) { $a = $Mk{$key_d_ex}[$idx] } }
  }
  if ( $a ) { $Mk{$key_ex}[$idx] =  $a }
# if not defined, set linker to first defined compiler
  unless ( $Mk{Linker}[$idx] ) { $Mk{Linker}[$idx] = $Mk{$key_ex}[$idx] }
 }

#-------------------------------------------------------------------------------
# define optional flags
#-------------------------------------------------------------------------------

 if ( $Mk{Source_Comment_Character}[$idx] ) {
  if ( $Mk{Source_Comment_Character}[$idx] eq "false" ) {
   undef $Mk{Source_Comment_Character}[$idx] 
  }
 }
 if ( $Mk{Run_File_Local}[$idx] ) {
  if ( $Mk{Run_File_Local}[$idx] eq "false" ) {
   undef $Mk{Run_File_Local}[$idx]
  }
 }
 unless ( $Mk{No_Warnings}[$idx] ) { $Mk{No_Warnings}[$idx] = $No_Warnings }
 unless ( $Mk{No_Warnings}[$idx] ) { $Mk{No_Warnings}[$idx] = $No_Warnings }
 unless ( $Mk{Auto_Update}[$idx] ) { $Mk{Auto_Update}[$idx] = $Auto_Update }
 unless ( $Mk{Preprocessor_Code}[$idx] ) {
  $Mk{Preprocessor_Code}[$idx] = $Preprocessor_Code
 }
 unless ( $Mk{Preprocessor_Compile}[$idx] ) {
  $Mk{Preprocessor_Compile}[$idx] = $Preprocessor_Compile
 }

}

#-------------------------------------------------------------------------------
# Reduce the size of Mk and add any update directories
#-------------------------------------------------------------------------------

undef %MkH;
$n = 0;
for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
 if ( -e $Mk{Source_Directory}[$idx] ) {
  $a = $Mk{Source_Directory}[$idx]; $a =~ s/$Initial_Directory//g;
  if ( -e $Mk{Updates_Level}[$idx].$a ) {
   foreach $var ( keys %Mk ) {
    $MkH{$var}[$n] = $Mk{$var}[$idx];
   }
   $MkH{Source_Directory}[$n] = $Mk{Updates_Level}[$idx].$a;
   $n++;
  }
  $a = $Mk{Source_Directory}[$idx]; $a =~ s/$Mk{Executable_Directory}[$idx]//g;
  if ( -e $Mk{Updates_Level}[$idx].$a ) {
   foreach $var ( keys %Mk ) {
    $MkH{$var}[$n] = $Mk{$var}[$idx];
   }
   $MkH{Source_Directory}[$n] = $Mk{Updates_Level}[$idx].$a;
   $n++;
  }
  $a = $Mk{Source_Directory}[$idx]; $a =~ s/$Mk{Version_Directory}[$idx]//g;
  if ( -e $Mk{Updates_Level}[$idx].$a ) {
   foreach $var ( keys %Mk ) {
    $MkH{$var}[$n] = $Mk{$var}[$idx];
   }
   $MkH{Source_Directory}[$n] = $Mk{Updates_Level}[$idx].$a;
   $n++;
  }
  foreach $var ( keys %Mk ) {
   $MkH{$var}[$n] = $Mk{$var}[$idx];
  }
  $n++;
 }
}
%Mk = %MkH;

#-------------------------------------------------------------------------------
# copy an input files
#-------------------------------------------------------------------------------

for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
 if ( -e $Mk{Executable_Directory}[$idx] && $Mk{Input_File}[$idx] ) {
  if ( ! -e $Mk{Input_File}[$idx] ) {
   $a = $Mk{Input_File}[$idx]; $a =~ s/.*\///g;
   $a = "$Home/.mk/$a";
   if ( -e "$a" ) { copy $a, $Mk{Input_File}[$idx] }
  }
  if ( ! -e $Mk{Input_File}[$idx] ) {
   $a = $Mk{Input_File}[$idx]; $a =~ s/.*\///g;
   $a = "$Mk{Version_Directory}[$idx]/$a";
   if ( -e "$a" ) { copy $a, $Mk{Input_File}[$idx] }
  }
  if ( ! -e $Mk{Input_File}[$idx] ) {
   $a = $Mk{Input_File}[$idx]; $a =~ s/.*\///g;
   $a = "$Mk{Version_Directory}[$idx]/run/$a";
   if ( -e "$a" ) { copy $a, $Mk{Input_File}[$idx] }
  }
 }
}

#-------------------------------------------------------------------------------
# if required, write help
#-------------------------------------------------------------------------------

 unless ( $Cmd ) { $Cmd = " " }
# write help to the screen
 if ( $Cmd eq " " ) { print $text1; exit }
 if ( $Cmd eq "h" || $Cmd eq "H" ) {
# if a help $String is present and not set to $Version, look for commands
  if ( $String ) {
   $Version =~ s/.*\///;
   if ( "$String" ne "$Version" ) {
#    print "$line";
#    print "mk help String: $String\n";
    print "$line";
    undef $a;
#   create a lower case version of $String
    $string = lc $String;
    @text2 = @extension_list = split(/\n/, $text2 );
#   search through each line of $text2
    foreach $text (@text2) {
     $_ = lc $text;
#    skip if $a is not defined yet
     if ( $a ) {
      if ( ! /^ / && /$string/ ) { $a = $string }
      elsif ( ! /^ / ) { $a = " " }
      if ( $a eq $string ) { print "$text\n" }
     }
#    define $a after finding a line starting with Alphabetical
     if ( /^alphabetical/ ) { $a = " " }
    }
    print "\nFor current settings, grep or look for the variable in the";
    print " file listed below.\n";
    print "$line";
   }
   else { print $text1 }
  }
  else { print $text1 }

  $list = " ";
# write help to the log file
  for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
   if ( index ( $list, $Mk{Log_File}[$idx] ) < 0 ) {
    print "More help text written to $Mk{Log_File}[$idx]\n";
    open OUT, ">$Mk{Log_File}[$idx]";
    print OUT $text1;
    print OUT $text2;
    @sortedkeys = sort ( keys(%Mk) );
#   write current settings
    print  OUT "Setup for Directory: $Initial_Directory\n";
    for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
     print OUT $line;
     print  OUT "  Source Directory Number: $idx\n";   
     foreach $var (@sortedkeys) {
      if ( $Mk{$var}[$idx] ) {
       write_out ( split( / /,"$var($idx): $Mk{$var}[$idx]" ) );
      }
     }
    }
#   write a list of all available model options
    print OUT $line;
    print OUT "Available options:\n";
    $n=0;
    for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
     $extension = "$Mk{Source_Extension}[$idx] $Mk{Include_Extension}[$idx]";
     @extension_list = split(/ /, $extension );
     foreach $extension (@extension_list) {
      foreach (<$Mk{Source_Directory}[$idx]/*.$extension>) {
       open IN, "$_";
       LINE: 
        while (<IN>) {
         next LINE if !/if/;
         s/\s+/ /g; s/\n/ /g; s/^ +//; next LINE if !/^#/;
         s/#/ /; next LINE if /^ include/i; next LINE if /^ endif/i;
         s/!/ /g; s/\(/ /g; ; s/\)/ /g; s/ ifndef/ /gi; s/ ifdef/ /gi;
         s/ elif/ /gi; s/ if/ /gi; s/ defined/ /gi; s/ define/ /gi;
         s/ endif/ /gi; s/\)//g; s/\(//g; s/\s+//g; s/\|\|/ /g; s/&&/ /g;
         push(@options, split(/\s+/, $_));
        }
       close IN;
      }
      $n++; 
     }
    }
    @options = mk_unique(sort(@options));
    write_out ( @options );
    $n = $#options + 1;
    print OUT "Number of available options: $n\n";
    close OUT;
   }
   $list = $list." ".$Mk{Log_File}[$idx]." ";
  }
  exit;
 }

#-------------------------------------------------------------------------------
# create code and data directories, if they do not exist
#-------------------------------------------------------------------------------

for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
 if ( $Mk{Code_Directory}[$idx] ) {
  unless ( -e $Mk{Code_Directory}[$idx] ) {
   mkdir ( $Mk{Code_Directory}[$idx],0777 );
   if ( $Mk{Uncompress_Command}[$idx] && $Mk{Uncompress_Extension}[$idx] ) {
    if ( -e "$Mk{Code_Directory}[$idx].$Mk{Compress_Extension}[$idx]" ) {
     chdir $Mk{Code_Directory}[$idx]; chdir "../";
     $a = $Mk{Code_Directory}[$idx]; $a =~ s/\/+$//; $a =~ s/.*\///g;
     $a = "$a.$Mk{Compress_Extension}[$idx]";
     system "$Mk{Uncompress_Command}[$idx] $a";
    }
   }
  }
 } 
 if ( $Mk{Backup_Directory}[$idx] ) {
  unless ( -e $Mk{Backup_Directory}[$idx] ) {
   mkdir ( $Mk{Backup_Directory}[$idx],0777 ) 
  }
 }
 if ( $Mk{Dependency_Directory}[$idx] ) {
  unless ( -e $Mk{Dependency_Directory}[$idx] ) {
   mkdir ( $Mk{Dependency_Directory}[$idx],0777 )
  }
 }
 if ( $Mk{Object_Directory}[$idx] ) {
  unless ( -e $Mk{Object_Directory}[$idx] ) {
   mkdir ( $Mk{Object_Directory}[$idx],0777 )
  }
 }
 if ( $Mk{Data_Directory}[$idx] ) {
  unless ( -e $Mk{Data_Directory}[$idx] ) {
   mkdir ( $Mk{Data_Directory}[$idx],0777 )
  }
 }
}

#-------------------------------------------------------------------------------
# clear the code directories, if required
#-------------------------------------------------------------------------------

if ( $Cmd eq "c" or $Cmd eq "C") {
 $list = " ";
 for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
  if ( index($list, $Mk{Code_Directory}[$idx]) < 0 ) {
   use File::Copy;
   if ( $Mk{Backup_Directory}[$idx] ) {
    @files = <$Mk{Backup_Directory}[$idx]/*>; if (@files) { unlink @files }
   }
   @files = <$Mk{Dependency_Directory}[$idx]/*>; if (@files) { unlink @files };
   @files = <$Mk{Object_Directory}[$idx]/*>; if (@files) { unlink @files };
   @files = <$Mk{Code_Directory}[$idx]/*>;
   if ( @files && $Mk{Backup_Directory}[$idx] ) {
    foreach $file (@files) { 
     if ($file ne $Mk{Backup_Directory}[$idx] ){
      copy ($file, $Mk{Backup_Directory}[$idx])
     }
    }
   }
   @files = <$Mk{Code_Directory}[$idx]/*>;if (@files) { unlink @files };
   print "Cleared $Mk{Code_Directory}[$idx]\n";
   $list = $list." ".$Mk{Code_Directory}[$idx]." ";
   if ( $Mk{Executable_File}[$idx] ) {
    if ( -e $Mk{Executable_File}[$idx] ) { unlink $Mk{Executable_File}[$idx] }
   }
  }
 }
 exit;
}

#-------------------------------------------------------------------------------
# check history file for any changes
#-------------------------------------------------------------------------------

# set Equality_Character to "<=" for History files
$Equality_Character = "<=";
$Changes = "false";

@sortedkeys = sort ( keys(%Mk) );
for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) { 
 if ( -e $Mk{Log_File}[$idx] ) { unlink $Mk{Log_File}[$idx] };
}

#-------------------------------------------------------------------------------
# create History file in each Dependency_Directory if it is missing
#-------------------------------------------------------------------------------

$list = " ";
for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {

 if ( ! $skip0 and index ( $list, $Mk{Dependency_Directory}[$idx] ) < 0 ) {

  $n = 0;
  $extension = "$Mk{Include_Extension}[$idx] $Mk{Source_Extension}[$idx]";
  @extension_list = split( / /, $extension);
  EXTENSION:
   foreach $extension (@extension_list) {
    if ( <$Mk{Code_Directory}[$idx]/*.$extension> ) { $n++; last EXTENSION}
   }
  

  if ( $Mk{No_Warnings}[$idx] eq "true" ) { open OUT, ">$Mk{Log_File}[$idx]" }
  else { open OUT, ">-" }
  if ( $n > 0 && ! -e "$Mk{Dependency_Directory}[$idx]/History" ) {
   $Changes = "true";
   print OUT "\n=>Warning: Code exists in $Mk{Code_Directory}[$idx] but ";
   print OUT "\n           the file $Mk{Dependency_Directory}[$idx]/History is";
   print OUT "\n           missing. If you continue, a new History file will";
   print OUT "\n           be written. This may make the code inconsistent.";
   print OUT "\n           This will apply to all source directories.\n";
   
   unless ( $Mk{No_Warnings}[$idx]  eq "true" ) {
    print OUT "\n  Continue (y/n)? "; $a = <STDIN>; chomp($a); $n = 0;
    if ( $a eq "n" ) { exit };
    $skip0 = "true";
   }
  }
  close OUT;
  
  if ( $n == 0 ) {
   open OUT, ">$Mk{Dependency_Directory}[$idx]/History";
   $Changes = "true";
   print OUT "    Equality_Character $Equality_Character\n";
   for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
    foreach $var (@sortedkeys) {
     if ( $Mk{$var}[$idx] ) {
      @out = split ( / /, "$var($idx) $Equality_Character $Mk{$var}[$idx]" );
      write_out ( @out );
     }
    }
   }
   close OUT;
  }
  
 }
 $list = $list." $Mk{Dependency_Directory}[$idx] ";
}


$list = " ";
for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {

 if ( index ( $list , $Mk{Code_Directory}[$idx] ) < 0 ) {
  undef %MkH; read_mk ( "$Mk{Dependency_Directory}[$idx]/History", %MkH );

#-------------------------------------------------------------------------------
# look for changes in source directories that use this code directory
#-------------------------------------------------------------------------------

  $n = 0;
  DIRECTORY:
   for $i ( 0 .. $#{ $Mk{Source_Directory} } ) {
    if ( $#{ $Mk{Source_Directory} } != $#{ $MkH{Source_Directory} } ) {
     $n++; last DIRECTORY
    }
    if ( $Mk{Code_Directory}[$idx] eq $MkH{Code_Directory}[$i] ) {
     if ( $Mk{Source_Directory}[$i] ne $MkH{Source_Directory}[$i] ) {
      $n++; last DIRECTORY
     }
    }
   }
   
  if ( $Mk{No_Warnings}[$idx] eq "true" ) { open OUT, ">>$Mk{Log_File}[$idx]" }
  else { open OUT, ">-" }

  if ( ! $skip1 and $n > 0 ) {  
   $Changes = "true";
   print OUT "\n=>Warning: Code Directory: $Mk{Code_Directory}[$idx]";
   print OUT "\n           Source Directories have changed.\n\n";
   print OUT "  Original:\n";
   for $i ( 0 .. $#{ $MkH{Source_Directory} } ) {
    print OUT "    $i: $MkH{Source_Directory}[$i]\n";
   }
   print OUT "  Current: \n";
   for $i ( 0 .. $#{ $Mk{Source_Directory} } ) {
    print OUT "    $i: $Mk{Source_Directory}[$i]\n" 
   }
  }
  
  if ( ! $skip1 and $n > 0 ) {  
   $a = "y";
   $Changes = "true";
   unless ( $Mk{No_Warnings}[$idx] eq "true" ) {
    print "\n  Remove ALL code in $Mk{Code_Directory}[$idx] (y/n)? ";
    $a = <STDIN>; chomp($a);
   }
   if ($a eq "y") {
    use File::Copy;
    if ( -e $Mk{Executable_File}[$idx] ) { unlink $Mk{Executable_File}[$idx] }
    if ( $Mk{Backup_Directory}[$idx] ) {
     @files = <$Mk{Backup_Directory}[$idx]/*>; if (@files) { unlink @files }
    }
    @files = <$Mk{Dependency_Directory}[$idx]/*>; if (@files) { unlink @files };
    @files = <$Mk{Object_Directory}[$idx]/*>; if (@files) { unlink @files };
    @files = <$Mk{Code_Directory}[$idx]/*>;
    if ( @files && $Mk{Backup_Directory}[$idx] ) {
     foreach $file (@files) { copy ( $file, $Mk{Backup_Directory}[$idx] ) }
    }
    @files = <$Mk{Code_Directory}[$idx]/*>; if (@files) { unlink @files };
    print OUT "Cleared $Mk{Code_Directory}[$idx]\n"; undef %MkH; 
    close OUT; 
    open OUT, ">$Mk{Dependency_Directory}[$idx]/History";
    $Changes = "true";
    for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
     foreach $var (@sortedkeys) {
      if ( $Mk{$var}[$idx] ) {
       @out = split( / /, "$var($idx) $Equality_Character $Mk{$var}[$idx]" );
       write_out ( @out );
      }
     }
    }
    close OUT;
   } else {
    print OUT "\n=>Warning: Affected code has not been removed. This may";
    print OUT "\n           make the code inconsistent. This will apply";
    print OUT "\n           to all source directories.\n";
    print OUT "\n  Continue (y/n)? "; $a = <STDIN>; chomp($a);
    if ( $a eq "n" ) { exit }; 
    $skip1 = "true"; 
    close OUT;
   }
  }
  else { close OUT }
 }
 $list = $list." $Mk{Code_Directory}[$idx] ";
}

for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {

 foreach $var ( keys(%MkH) ) {

  if ( $Mk{$var}[$idx] ) {
   $val = $Mk{$var}[$idx]; $val =~ s/ +$// ; $val =~ s/^ +//;
  }
  else { $val = " " }
  if ( $MkH{$var}[$idx] ) {
   $valH = $MkH{$var}[$idx]; $valH =~ s/ +$// ; $valH =~ s/^ +//;
  }
  else { $valH = " " }
  
  if ( $val ne $valH ) {
   $Changes = "true";
  
#-------------------------------------------------------------------------------
#  look for model option changes
#-------------------------------------------------------------------------------

   if ( $var eq "Model_Options" ) {
    $option_changes = "";
    
    unless ($MkH{$var}[$idx]) {$MkH{$var}[$idx] = " "}
    unless ($Mk{$var}[$idx]) {$Mk{$var}[$idx] = " "}
    @new_options = split(/ +/, "$Mk{$var}[$idx]");
    foreach (@new_options) {
     if (index(" $MkH{$var}[$idx] ", " $_ ") < 0) {
      $option_changes = "$option_changes $_";
     }
    }
    @old_options = split(/ +/, "$MkH{$var}[$idx]");
    foreach (@old_options) {
     if (index(" $Mk{$var}[$idx] ", " $_ ") < 0) {
      $option_changes = "$option_changes $_";
     }
    }
    $option_changes =~ s/^ +//; $option_changes =~ s/ +$//;
    if ($option_changes) {
     undef @file_list;
     @option_changes = split(/ +/, $option_changes);
     foreach $option (@option_changes) {
      $n = 0;
      foreach $file (<$Mk{Dependency_Directory}[$idx]/*.d>) {
       open IN, "$file";
       while (<IN>) {
        chomp($_);
	if ( index (" $_ ", " $option ") > 0 ) {
         $file =~ s/.*\///; $file =~ s/\.d$//;
	 if ( -e "$Mk{Code_Directory}[$idx]/$file" ) { push(@file_list, $file) }
	}
       }
      }
     }
     if ( ! $skip2 and @file_list ) {
      @file_list = mk_unique(sort(@file_list));     
      if ( $Mk{No_Warnings}[$idx] eq "true" ) {
       open OUT, ">>$Mk{Log_File}[$idx]";
      }
      else { open OUT, ">-" }
      print OUT "\n=>Warning: Source_Directory: $Mk{Source_Directory}[$idx]\n";
      print OUT "\n  The following options have been added or removed:\n";
      write_out ( split( / /, "$option_changes" ) );      
      print OUT "\n  The following files are affected:\n";
      write_out ( @file_list );      
      $a = "y";
      unless ( $Mk{No_Warnings}[$idx] eq "true" ) {
       print "\n  Remove affected code (y/n)? "; $a = <STDIN>; chomp($a);
      }
      if ($a eq "y") {
       print OUT "\n  Deleting\n";
       write_out ( @file_list );
       chdir $Mk{Code_Directory}[$idx];
       if ( -e $Mk{Executable_File}[$idx] ) {
        unlink $Mk{Executable_File}[$idx]
       }
       unlink @file_list;
      }
      else {
       print OUT "\n=>Warning: Affected code has not been removed. Any newly";
       print OUT "\n           made code will use the new options. This may";
       print OUT "\n           make the code inconsistent. This will apply";
       print OUT "\n           to all source directories.\n";
       print OUT "\n  Continue (y/n)? "; $a = <STDIN>; chomp($a);
       if ( $a eq "n" ) { exit };
       $skip2 = "true"; 
      }
      close OUT;
     }
    }
    $var = "";
   }

#-------------------------------------------------------------------------------
#  look for compiler changes
#-------------------------------------------------------------------------------

   if ( $var ) {
    @extension_list = split(/ /,$Mk{Source_Extension}[$idx]);
    foreach $extension (@extension_list) {
     $key = "Compiler_".$extension;
     if ( $var eq $key ) {
      @files = <$Mk{Object_Directory}[$idx]/*>; 
      if (! $skip3 and @files) {
       if ( $Mk{No_Warnings}[$idx] eq "true" ) {
        open OUT, ">>$Mk{Log_File}[$idx]";
       }
       else { open OUT, ">-" }
       print OUT "\n=>Warning: Source_Directory: $Mk{Source_Directory}[$idx]";
       print OUT "\n           $key for this source directory has changed.\n";
       print OUT "\n  Original $var: $valH ";
       print OUT "\n  Current $var:  $val \n";
       $a = "y";
       unless ( $Mk{No_Warnings}[$idx] eq "true" ) {
        print OUT "\n  Remove code from $Mk{Object_Directory}[$idx] (y/n)? ";
        $a = <STDIN>; chomp($a);
       }
       if ($a eq "y") { 
        if ( -e $Mk{Executable_File}[$idx] ) {
         unlink $Mk{Executable_File}[$idx]
        }
        print OUT "\n  Deleting\n";
        write_out ( @files );
        unlink @files;
       }
       else {
	print OUT "\n=>Warning: Affected code has not been removed. This may";
	print OUT "\n           make the code inconsistent. This will apply";
        print OUT "\n           to all source directories.\n";
	print OUT "\n  Continue (y/n)? "; $a = <STDIN>; chomp($a);
	if ( $a eq "n" ) { exit };
        $skip3 = "true"; 
       }
       close OUT;
      }
      $var = "";
     }
    }
   }

#-------------------------------------------------------------------------------
#  look for linker or library changes
#-------------------------------------------------------------------------------

   if ( $var && $Mk{Executable_File}[$idx] ) {
    if ( $var eq Linker || $var eq Libraries  || $var eq Libraries_Directory  ||
         $var eq Libraries_Directory_Prefix ) {
     if ( ! $skip4 and -e $Mk{Executable_File}[$idx] ) {
      if ( $Mk{No_Warnings}[$idx] eq "true" ) {
       open OUT, ">>$Mk{Log_File}[$idx]";
      }
      else { open OUT, ">-" }
      print OUT "\n=>Warning: Source_Directory: $Mk{Source_Directory}[$idx]";
      print OUT "\n           Linking for this source directory has changed.\n";
      print OUT "\n  Original $var: $valH";
      print OUT "\n  Current $var:  $val\n";
      $a = "y";
      unless ( $Mk{No_Warnings}[$idx] eq "true" ) {
       print OUT "\n  Remove $Mk{Executable_File}[$idx] (y/n)? ";
       $a = <STDIN>; chomp($a);
      }
      if ($a eq "y") { unlink $Mk{Executable_File}[$idx] }
      else {
       print OUT "\n=>Warning: Affected code has not been removed. This may";
       print OUT "\n           make the code inconsistent. This will apply";
       print OUT "\n           to all source directories \n";
       print OUT "\n  Continue (y/n)? "; $a = <STDIN>; chomp($a);
       if ( $a eq "n" ) { exit };
       $skip4 = "true";
      }
      close OUT;
     }
     $var = ""
    }
   }

#-------------------------------------------------------------------------------
#  warn of other changes
#-------------------------------------------------------------------------------

   if ( $var ) {
    open OUT, ">>$Mk{Log_File}[$idx]";
    print OUT "\n=>Warning: Source_Directory: $Mk{Source_Directory}[$idx]\n";
    print OUT "\n  Original $var: $valH ";
    print OUT "\n  Current $var:  $val \n";
   }

  }

 }
 
#-------------------------------------------------------------------------------
# auto update
#-------------------------------------------------------------------------------

 if ( $Mk{Auto_Update}[$idx] eq "true" ) {
  $extension = "$Mk{Include_Extension}[$idx] $Mk{Source_Extension}[$idx]";
  @extension_list = split( / /, $extension );
  foreach $extension (@extension_list) {
   @file_list = <$Mk{Source_Directory}[$idx]/*.$extension>;
   foreach $file (@file_list) {
    $time_s = 0; if ( -e $file) { $time_s = (stat($file))[9] }
    $file_s = $file;
    $file =~ s/.*\///g; $file = "$Mk{Code_Directory}[$idx]/$file";
    $time_c = 0; if ( -e $file) { $time_c = (stat($file))[9] }
    if ( ! -e $file) {
     $ext = $file; $ext =~ s/.*\.//g; $ext =~ s/ //g;
     if ( index( $Mk{Source_Extension}[$idx], $ext ) > -1 ) {
      $file =~ s/.*\///g; $file =~ s/\..*//;
      $file = "$Mk{Code_Directory}[$idx]/$file.$Mk{Code_Extension}[$idx]";
      $time_c = 0; if ( -e $file) { $time_c = (stat($file))[9] } 
     }
    }
    if ( -e $file && $time_c lt $time_s ) {
     $Changes = "true";
     if ( $Mk{No_Warnings}[$idx] eq "true" ) {
      open OUT, ">>$Mk{Log_File}[$idx]";
     }
     else { open OUT, ">-" }
     print OUT "\n=>Warning: Source_Directory: $Mk{Source_Directory}[$idx]\n";
     print OUT "           $file is older than\n";
     print OUT "           $file_s\n";
     $a = "y";
     unless ( $Mk{No_Warnings}[$idx] eq "true" ) {
      print "\n  Remove $file (y/n)? "; $a = <STDIN>; chomp($a);
     }
     if ($a eq "y") {
      print OUT "  Deleting $file\n";
      unlink $file;
      if ( -e $Mk{Executable_File}[$idx] ) {
       unlink $Mk{Executable_File}[$idx]
      }
     }
    }  
   }
  }
 }
 
#-------------------------------------------------------------------------------
# check Preprocessor_Compile
#-------------------------------------------------------------------------------

 if ( $Mk{Preprocessor_Compile}[$idx] eq "true" && $Mk{Preprocessor}[$idx] ) {
  @extension_list = split( / /, $Mk{Source_Extension}[$idx] );
  foreach $extension (@extension_list) {
   $key = "Compiler_$extension";
   if ( $Mk{$key}[$idx] && $Mk{Preprocessor}[$idx] ) {
    if ( index( $Mk{$key}[$idx], $Mk{Preprocessor}[$idx] ) < 0 ) {
     $Changes = "true";
#    if Preprocessor not found in Compiler undefine Preprocessor
     open OUT, ">>$Mk{Log_File}[$idx]";
     print OUT "\n=>Warning: Preprocessor_Compile is set to true \n";
     print OUT "           but $Mk{Preprocessor}[$idx] not found in \n";
     print OUT "           $Mk{$key}[$idx].\n";
     print OUT "           Undefining Preprocessor\n";
     close OUT;
     undef $Mk{Preprocessor}[$idx];
    }
   }
  }
 }
 
}

#-------------------------------------------------------------------------------
# write History file
#-------------------------------------------------------------------------------

$list = " ";
for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) { 
 if ( index ( $list, $Mk{Dependency_Directory}[$idx] ) < 0 ) {
  open OUT, ">$Mk{Dependency_Directory}[$idx]/History";
  print OUT "    Equality_Character $Equality_Character\n";
  for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
   foreach $var (@sortedkeys) {
    if ( $Mk{$var}[$idx] ) {
     @out = split ( / /, "$var($idx) $Equality_Character $Mk{$var}[$idx]" );
     write_out ( @out );
    }
   }
  }
  close OUT;
 }
}

#-------------------------------------------------------------------------------
# make mk_include_file in each code directory
#-------------------------------------------------------------------------------

$list = " ";
for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
 if ( $Mk{Mk_Include_File}[$idx] ) {
  $a = "$Mk{Code_Directory}[$idx]/$Mk{Mk_Include_File}[$idx]";
  if ( ! -e $a ) { $Changes = "true" }
 }
 if ( $Mk{Mk_Include_File}[$idx] && $Changes eq "true" ) {
  if ( index( $list, $Mk{Code_Directory}[$idx] ) < 0 ) {
   chdir $Mk{Code_Directory}[$idx];
   open OUT, ">$Mk{Mk_Include_File}[$idx]";
   for $i ( 0 .. $#{ $Mk{Source_Directory} } ) {
    foreach $var (@sortedkeys) {
     if ( $Mk{$var}[$i] ) {
      @out = split ( / /, "$var($i) $Equality_Character $Mk{$var}[$i]" );
      $c = 0; $s = "	"; $w = 80;
      if ( $Mk{Mk_Include_Start}[$idx] ) { $s = $Mk{Mk_Include_Start}[$idx] }
      if ( $Mk{Mk_Include_Width}[$idx] ) { $w = $Mk{Mk_Include_Width}[$idx] }
      print OUT "$s\'"; $c += length($s);
      foreach (@out) {
       s/^ +//; s/$ +//;
       $c += length(" $_");
       if ($c > $w) { print OUT "\'\n$s\'"; $c = length("$s\' $_") }
       if ($c <= $w) { print OUT " $_" }
       else { print OUT " \'\n$s\' variable is too long to print" }
      }
      print OUT "\'\n";
     }
    }
   }
   close OUT;
   $list = "$list $Mk{Code_Directory}[$idx] ";
  }
 }
}


#-------------------------------------------------------------------------------
# add search path to libraries
#-------------------------------------------------------------------------------

for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
 if ( $Mk{Libraries_Directory}[$idx] ) {
  undef $a;
  @list = split ( /\s+/, $Mk{Libraries_Directory}[$idx] );
  foreach (@list) {
   if ($a) {$a = "$a $Mk{Libraries_Directory_Prefix}[$idx]$_"}
   else {$a = "$Mk{Libraries_Directory_Prefix}[$idx]$_"}
  }
  if ($a && $Mk{Libraries}[$idx]) {
   $Mk{Libraries}[$idx] = "$a $Mk{Libraries}[$idx]" 
  }
 }
}


#-------------------------------------------------------------------------------
# add Option_Prefix to Model_Options
#-------------------------------------------------------------------------------

for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
 undef @model_options;
 @model_options = split ( /\s+/, $Mk{Model_Options}[$idx] );
 $Mk{Model_Options}[$idx] = "";
 foreach $option (@model_options) {
  $option = "$Mk{Option_Prefix}[$idx]$option";
  $Mk{Model_Options}[$idx] = "$Mk{Model_Options}[$idx] $option ";
 }
# if using Preprocessor_Compile, find Preprocessor in Compiler and add options
 if ( $Mk{Preprocessor_Compile}[$idx] eq "true" &&  $Mk{Preprocessor}[$idx] ) {
  $extension = "$Mk{Source_Extension}[$idx]";
  @extension_list = split( / /, $extension);
  foreach $extension (@extension_list) {
   $key = "Compiler_$extension";
   if ( $Mk{$key}[$idx] ) {
    $a = $Mk{Preprocessor}[$idx];
    $Mk{$key}[$idx] =~ s/$a/$a $Mk{Model_Options}[$idx]/;
   }
  }
 }
}

#-------------------------------------------------------------------------------
# make code
#-------------------------------------------------------------------------------

if ( $Cmd eq "f" ) { mk_all_code () }
if ( $Cmd eq "o" ) { $Cmd = "f"; mk_all_code (); $Cmd = "o"; mk_all_code () }

#-------------------------------------------------------------------------------
# make executable
#-------------------------------------------------------------------------------

elsif ( $Cmd eq "e" || $Cmd eq "r" ||  $Cmd eq "s" || $Cmd eq "q" ) {

 $Cmd = "f"; mk_all_code (); $Cmd = "o"; mk_all_code (); $Cmd = $ARGV[0];
 
 $list_exe = " ";
 for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {

  chdir $Mk{Executable_Directory}[$idx];

  if ( $Mk{Executable_File}[$idx] && $Mk{Linker}[$idx] ) {
   if ( ! -e $Mk{Executable_File}[$idx]) {
    $list = " ";
    for $i ( 0 .. $#{ $Mk{Source_Directory} } ) {
     if ( $Mk{Executable_File}[$idx] eq $Mk{Executable_File}[$i] ) {
      if ( index( $list, $Mk{Object_Directory}[$idx] ) < 0 ) { 
       $list = "$list $Mk{Object_Directory}[$i]/*.$Mk{Object_Extension}[$i] ";
      }
     }
    }
    if ( index ( $list_exe, $Mk{Executable_File}[$idx] ) < 0 ) { 
     if ( $Mk{Libraries}[$idx] ) { $list = "$list $Mk{Libraries}[$idx]" }
     $a = "$Mk{Linker}[$idx] $Mk{Executable_File}[$idx] $list "; $a =~ s/ +/ /g;
     system "$a >> $Mk{Log_File}[$idx] 2>&1";
     open OUT, ">>$Mk{Log_File}[$idx]";
     print OUT "Executable Directory: $Mk{Executable_Directory}[$idx]\n";
     close OUT;
     open OUT, ">>$Mk{Log_File}[$idx]";
     $a =~ s/$Mk{Executable_Directory}[$idx]\///g;
     print OUT "Linking: $a\n";
     close OUT;
     $list_exe = "$list_exe $Mk{Executable_File}[$idx] ";
    }
    if ( -e $Mk{Executable_File}[$idx] ) {
     open OUT, ">>$Mk{Log_File}[$idx]";
      print OUT "Made $Mk{Executable_File}[$idx]\n";
      close OUT 
     }
   }
  }
  else { 
   open OUT, ">>$Mk{Log_File}[$idx]";
   print OUT "\nSource_Directory: $Mk{Source_Directory}[$idx]\n";
   print OUT "Executable_File or Linker is not defined.\n";
   print OUT "Executable file was not created.\n";
   close OUT 
  }

#-------------------------------------------------------------------------------
# copy over data
#-------------------------------------------------------------------------------

  $list = " ";
  if ( $Mk{Data_Directory}[$idx] && $Mk{Data_Source}[$idx] ) {
   if ( index( $list, $Mk{Data_Directory}[$idx] ) < 0 ) { 

#   copy any updated data first
    $a = $Mk{Data_Source}[$idx]; $a =~ s/.*$Version\///g;
    $a = "$Mk{Updates_Level}[$idx]/$a";
    if ( -e $a ) {
     @files = <$a/*>;
     foreach $file (@files) {
      if ( -e $file) {
       $a = $file; $a =~ s/.*\///g;
       if ( ! -e "$Mk{Data_Directory}[$idx]/$a" ){
        copy $file, "$Mk{Data_Directory}[$idx]/$a";
       }
      }
     }
    }

    if ( -e $Mk{Data_Source}[$idx] ) {
     @files = <$Mk{Data_Source}[$idx]/*>;
     foreach $file (@files) { 
      if ( -e $file) {
       $a = $file; $a =~ s/.*\///g;
       if ( ! -e "$Mk{Data_Directory}[$idx]/$a" ){
        copy $file, "$Mk{Data_Directory}[$idx]/$a";
       }
      }
     }
    }
    
    $list = "$list $Mk{Data_Directory}[$idx] ";
   }
  }

 }

#-------------------------------------------------------------------------------
# look if a local machine Run_File_Local is set
#-------------------------------------------------------------------------------

 if ( $Cmd eq "r" ) {
  for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
   if ( $Mk{Run_File_Local}[$idx] ) {
    $Cmd = "q";
    $Mk{Run_File}[$idx] = $Mk{Run_File_Local}[$idx];
   }
  }
 }

#-------------------------------------------------------------------------------
# run on the local machine
#-------------------------------------------------------------------------------

 if ( $Cmd eq "r" ) {
  $list_exe = " ";
  for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
   chdir $Mk{Executable_Directory}[$idx];
   if ( $Mk{Executable_File}[$idx] ) {
    if ( -e $Mk{Executable_File}[$idx] ) {
     if ( index ( $list_exe, $Mk{Executable_File}[$idx] ) < 0 ) { 
      chdir $Mk{Executable_Directory}[$idx];
      $a = "touch $Mk{Executable_File}[$idx].run;";
      if  ( $Mk{Nice}[$idx] ) { $a = "$a nice $Mk{Nice}[$idx]" };
      $a = "$a $Mk{Executable_File}[$idx]";
      if ( $Mk{Output_File}[$idx] ) {
       unlink "$Mk{Output_File}[$idx]";
       $a = "$a > $Mk{Output_File}[$idx]";
      }
      $a = "$a 2>&1; rm $Mk{Executable_File}[$idx].run";
      system "( $a )&";
      $list_exe = $list_exe." $Mk{Executable_File}[$idx] ";
     }
    }
   }
  }
 }

#-------------------------------------------------------------------------------
# copy and run submit file 
#-------------------------------------------------------------------------------

 elsif ( $Cmd eq "s" || $Cmd eq "q" ) {
  $list_exe = " ";
  for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
   chdir $Mk{Executable_Directory}[$idx];
   if ( $Mk{Run_File}[$idx] ) {
    if ( -e $Mk{Executable_File}[$idx] ) {
     if ( index ( $list_exe, $Mk{Executable_File}[$idx] ) < 0 ) {
      $a = $Mk{Run_File}[$idx];
      unless ( -e $a ) {
       $a = "$Mk{Version_Directory}[$idx]/$a"
      }
      if ( -e $a ) {
       undef @out;
       open IN, $a; while (<IN>) { push ( @out, "$_" ) }; close IN;
       $a = "$Mk{Executable_File}[$idx].run";
       system "which perl > $a";
       open IN, $a; while (<IN>) { chomp($_); $perl = $_ }; close IN;
       open OUT, ">$a";
       print OUT "#! $perl -w\n";
       print OUT "\$Executable_Directory = \"$Mk{Executable_Directory}[$idx]\";\n";
       print OUT "\$Executable_File = \"$Mk{Executable_File}[$idx]\";\n";
       if ( $Mk{Number_Processors}[$idx] ) {
        print OUT "\$Number_Processors = \"$Mk{Number_Processors}[$idx]\";\n";
       }
       if ( $Mk{Input_File}[$idx] ) {
        print OUT "\$Input_File = \"$Mk{Input_File}[$idx]\";\n";
       }
       if ( $Mk{Output_File}[$idx] ) {
        print OUT "\$Output_File = \"$Mk{Output_File}[$idx]\";\n";
       }
       if ( $Mk{Step_Start_File}[$idx] ) {
        print OUT "\$Step_Start_File = \"$Mk{Step_Start_File}[$idx]\";\n";
       }
       if ( $Mk{Step_End_File}[$idx] ) {
        print OUT "\$Step_End_File = \"$Mk{Step_End_File}[$idx]\";\n";
       }
       if ( $Machine_Name ) {
        print OUT "\$Machine_Name = \"$Machine_Name\";\n";
       }
       if ( $Mk{Change_Mount}[$idx] ) {
        print OUT "\$Change_Mount = \"$Mk{Change_Mount}[$idx]\";\n";
       }
       if ( $Mk{Run_Directory}[$idx] ) {
        print OUT "\$Run_Directory = \"$Mk{Run_Directory}[$idx]\";\n";
       }
       if ( $Mk{Run_Copy}[$idx] ) {
        print OUT "\$Run_Copy = \"$Mk{Run_Copy}[$idx]\";\n";
       }
       if ( $Mk{Notify_User}[$idx] ) {
        if ( index ( $Mk{Notify_User}[$idx], "@" ) > -1 ) { 
         $Mk{Notify_User}[$idx] =~ s/\@/\\@/g;
        }  
        print OUT "\$Notify_User = \"$Mk{Notify_User}[$idx]\";\n";
       }
       if ( $Mk{End_String}[$idx] ) {
        print OUT "\$End_String = \"$Mk{End_String}[$idx]\";\n";
       }
       if ( $Mk{Nice}[$idx] ) {
        print OUT "\$Nice = \"$Mk{Nice}[$idx]\";\n";
       }
       foreach (@out) { print OUT "$_" };
       close OUT;
       system "chmod +x $a";
       system "$a $Cmd";
      }
      $list_exe = $list_exe." $Mk{Executable_File}[$idx] ";
     }
    }
   }
  }
 }
 
}

#-------------------------------------------------------------------------------
# make individual files
#-------------------------------------------------------------------------------

else {
 for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
  chdir $Mk{Code_Directory}[$idx];
  mk_code ( $Cmd );
 }
}

#-------------------------------------------------------------------------------
# compress code files if required and finish
#-------------------------------------------------------------------------------

$list = " ";
for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
 if ( index ( $list, $Mk{Code_Directory}[$idx] ) < 0 ) {
  if ( $Mk{Compress_Command}[$idx] && $Mk{Compress_Extension}[$idx] ) {
   if ( -e "$Mk{Code_Directory}[$idx]" ) {
    chdir "$Mk{Code_Directory}[$idx]"; chdir "../";
    $a = $Mk{Code_Directory}[$idx]; $a =~ s/\/+$//; $a =~ s/.*\///g;
    system "$Mk{Compress_Command}[$idx] $a.$Mk{Compress_Extension}[$idx] $a";
    if ( -e "$a.$Mk{Compress_Extension}[$idx]" ) {
     system "rm -rf $Mk{Code_Directory}[$idx]"
    }
   }
  }
 }
 $list = " $list $Mk{Code_Directory}[$idx] ";
}

$list = " ";
for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
 if ( index ( $list, $Mk{Log_File}[$idx] ) < 0 ) { 
  open OUT, ">>$Mk{Log_File}[$idx]"; print OUT "Done: mk $Cmd\n"; close OUT;
 }
 $list = " $list $Mk{Log_File}[$idx] ";
}
exit;

#*******************************************************************************
# subroutine mk_all_code
#*******************************************************************************
sub mk_all_code {

 my $file; my $extension; my $extension_list; my @file_list;

 for $idx ( 0 .. $#{ $Mk{Source_Directory} } ) {
  chdir $Mk{Code_Directory}[$idx];
  $extension = "$Mk{Include_Extension}[$idx] $Mk{Source_Extension}[$idx]";
  @extension_list = split( / /, $extension );
  foreach $extension (@extension_list) {
   @file_list = <$Mk{Source_Directory}[$idx]/*.$extension>;
   foreach $file (@file_list) {
    $file =~ s/.*\///g;
    mk_code ( $file );
   }
  }  
 }
 1;
}

#*******************************************************************************
# subroutine mk_code
#*******************************************************************************
sub mk_code {

 my $file; $file = $_[0]; $file=~ s/.*\///g;
 my $base; my $ext; my @dep; my $i; my $file_s;
 my $file_d; my $file_c; my $file_o; my @file_list; my $tmp; 
 my @extension_list; my $extension; my $n;

 $base = $file; $base =~ s/.*\///g; $base =~ s/\..*//;
 $ext = $file;  $ext =~ s/.*\///g; $ext =~ s/.*\.//;
 $tmp = "temporary_file";

# if file has a code or object extension find the source file
 DIRECTORY:
  for $i ( 0 .. $#{ $Mk{Source_Directory} } ) {
   if ( $Mk{Code_Extension}[$i] && $Mk{Object_Extension}[$i] ) {
    if ( $ext eq $Mk{Code_Extension}[$i] || $ext eq $Mk{Object_Extension}[$i] ) {
     if ( $Mk{Code_Directory}[$i] eq $Mk{Code_Directory}[$idx] ) {
      @extension_list = split( / /, $Mk{Source_Extension}[$idx] );
      foreach $extension (@extension_list) {
       if ( -e "$Mk{Source_Directory}[$idx]/$base.$extension" ) { 
        $file = "$base.$extension";
        $ext = $extension;
        last DIRECTORY;
       }
      }
     }
    }
   }
  }
 
 unless ( $made_list ) { $made_list = " " }
 if ( index($made_list, $Mk{Code_Directory}[$idx]."/".$file) > -1 ) { return } 
 
#-------------------------------------------------------------------------------
# create the code file
#-------------------------------------------------------------------------------
 
 $i = 0;
 DIRECTORY:
  for $n ( 0 .. $#{ $Mk{Source_Directory} } ) {

   $file_s = "$Mk{Source_Directory}[$i]/$file";
   $file_d = "$Mk{Dependency_Directory}[$i]/$file.d";
   $file_c = "$Mk{Code_Directory}[$i]/$file";
   if ( $Mk{Preprocessor}[$i] ) {
    if ( $Mk{Preprocessor_Code}[$i] ne "true" ) {
     if ($Mk{Preprocessor_Compile}[$i] ne "true" ) {
      if ( index($Mk{Source_Extension}[$i], $ext) > -1 ) {
       $file_c = "$Mk{Code_Directory}[$i]/$base.$Mk{Code_Extension}[$i]";
       $file_d = "$Mk{Dependency_Directory}[$i]/$base.$Mk{Code_Extension}[$i].d";
      }
     }
    }
   }
   
   if ( -e $file_s && $Mk{Code_Directory}[$i] eq $Mk{Code_Directory}[$idx] ) {
    if ( ! -e $file_c ) {
     if ( ! $Mk{Preprocessor}[$i] ) { copy $file_s, $file_c }
     if ( $Mk{Preprocessor_Code}[$i] eq "true" ) { copy $file_s, $file_c }
     if ( $Mk{Preprocessor_Compile}[$i] eq "true" ) { copy $file_s, $file_c }
#     if ( $ext eq $Mk{Code_Extension}[$i] ) {copy $file_s, $file_c }
     if ( ! -e $file_c ) {
#     copy the source file to $tmp and preprocess
      unlink $file_d;
      mk_dep ( $i, $file_s, $file_d );
      undef @out; read_cpp ( $i, $file_d );
      if ( $Mk{Source_Comment_Character}[$i] ) {
       push ( @out, "$Mk{Source_Comment_Character}[$i] source file: $file_s\n" )
      }
      open IN, $file_s; while (<IN>) { push ( @out, "$_" ) }; close IN;
      open OUT, ">$tmp"; foreach (@out) { print OUT "$_" }; close OUT;
#     remove $output and change preprocessor includes to code includes in $input
      unlink $file_c; undef @out;
      open IN, "$tmp";
      while (<IN>) {
       if (/include/) { s/</\"/; s/>/\"/ }
       if (/#/) {
        s/#.*include/      include/g; s/endif.*/endif/g; s/else.*/else/g
       }
       push (@out, $_);
      }
      close IN;
      open OUT, ">$tmp"; foreach (@out) { print OUT "$_" }; close OUT;
#     preprocess $tmp to $file_c
      $a = "$Mk{Preprocessor}[$i] $Mk{Model_Options}[$i]";
      system "$a < $tmp > $file_c 2>>  $Mk{Log_File}[$i]";
#     remove $tmp and replace multiple blank lines with single ones in $file_c
      unlink $tmp; undef @out; $a = "";
      open IN, "$file_c";
      while (<IN>) {
       chomp($_); s/^#.*//; s/\s+$//;
       if ( $_ || $a ne $_ ) { push ( @out, $_ ) }
       $a = $_;
      }
      close IN;
      open OUT, ">$file_c"; foreach (@out) { print OUT "$_\n" }; close OUT;
      unlink $tmp;
      unlink $file_d;
     }
     if ( -e $file_c ) { 
      open OUT, ">>$Mk{Log_File}[$idx]";
      print OUT "Made $file_c\n";
      close OUT
     }
    }
    if (! -e $file_d ) {
     mk_dep ( $i, $file_c, $file_d );
     if ( -e $file_d ) {
      open OUT, ">>$Mk{Log_File}[$idx]";
      print OUT "Made $file_d\n";
      close OUT
     } 
    }
    last DIRECTORY;  
   }
   $i++;

  }

 if (  $i > $#{ $Mk{Source_Directory} } ) { return }
 if ( ! -e $file_c ) {
  open OUT, ">>$Mk{Log_File}[$idx]";
  print OUT "\n=>Warning: $file_c could not be made.\n";
  close OUT;
  return;
 }

 if ( $Cmd eq $file_c || "$Mk{Code_Directory}[$i]/$Cmd" eq $file_c ) { return };
 if ( $Cmd eq "f" || $Cmd eq "F" ) { return };
 $a = $Cmd; $a =~ s/.*\.//g;
 if ( $a eq "f" || $a eq "F" ) { return };

#-------------------------------------------------------------------------------
# create the object file
#-------------------------------------------------------------------------------

 if ( index($Mk{Source_Extension}[$i], $ext) > -1 ) {

  $file_o = "$Mk{Object_Directory}[$i]/$base.$Mk{Object_Extension}[$i]";
  $time_o = 0; if ( -e $file_o) { $time_o = (stat($file_o))[9] }
  undef %Checked; 
  if ( -e $file_o ) {
   if ( check_time ( $i, $file_d, $time_o ) ) { unlink $file_o }
  }
  
  if (! -e $file_o ) {
  
   if ( $Mk{Executable_File}[$idx] ) { 
    if ( -e $Mk{Executable_File}[$idx] ) { unlink $Mk{Executable_File}[$idx] }
   }
   mk_dep ( $i, $file_c, $file_d );

   unless ( $Mk{Preprocessor_Compile}[$i] eq "true" ) {
    if ( $Mk{Preprocessor_Code}[$i] eq "true" && $Mk{Preprocessor}[$i] ) {
     if ( $ext ne $Mk{Code_Extension}[$i] ) {
      $tmp = $file_c;
      $file_c = "$Mk{Code_Directory}[$i]/$base.$Mk{Code_Extension}[$i]";
      $a = "$Mk{Preprocessor}[$i] $Mk{Model_Options}[$i] < $tmp > $file_c";
      system "$a 2>> $Mk{Log_File}[$idx]";
     }
    }
   }
   
   if ( $Mk{Module_Extension}[$i] ) {
    @file_list = <$Mk{Object_Directory}[$i]/*.$Mk{Module_Extension}[$i]>;
    foreach $f (@file_list) { copy $f, $Mk{Code_Directory}[$i] }
   }

   $key = "Compiler_$ext";
   if ( $Mk{$key}[$i] ) {
    system "$Mk{$key}[$i] $file_c >> $Mk{Log_File}[$i] 2>&1";
   }

   if ( $Mk{Module_Extension}[$i] ) {
    @file_list = <$Mk{Code_Directory}[$i]/*.$Mk{Module_Extension}[$i]>;
    foreach $f (@file_list) { copy $f, $Mk{Object_Directory}[$i]; unlink $f };
   }

   if ( $Mk{Object_Extension}[$i] ) {
    @file_list = <$Mk{Code_Directory}[$i]/*.$Mk{Object_Extension}[$i]>;
    foreach $f (@file_list) { copy $f, $Mk{Object_Directory}[$i]; unlink $f };
   }

   if ( $Mk{Preprocessor_Code}[$i] eq "true" && $Mk{Preprocessor}[$i] ) {
    if ( $ext ne $Mk{Code_Extension}[$i] ) {
     unlink $file_c;
    }
   }

   if ( -e $file_o ) {
    open OUT, ">>$Mk{Log_File}[$idx]";
    print OUT "Made $file_o\n";
    close OUT 
   }

  }
 }
 $a = $Mk{Code_Directory}[$i]."/".$file;
 $made_list = $made_list." $a ";
 1;
}


#*******************************************************************************
# subroutine check_time
#*******************************************************************************
sub check_time {

 my $i; $i = $_[0];
 my $file_d; $file_d = $_[1];
 my $time_o; $time_o = $_[2];
 my $out; my $time_d; my $file; my $n; my @in; 
 
 open IN, $file_d; 
 LINE:
  while (<IN>) {
   chomp $_; s/ +//g;
   last LINE if /#/;   
   push ( @in, $_ )
  }
 close IN;

 LINE:
  foreach (@in) {
   unless ( $Checked{$_} ) {
    $Checked{$_} =  $_;
    $Checked{$_} = 0; if ( -e $_ ) { $Checked{$_} = (stat($_))[9] }
    if ( ! -e $_ || $time_o < $Checked{$_} ) { $out = 1; last LINE}
    $file = "$Mk{Dependency_Directory}[$i]/$_.d";
    unless ( $Checked{"$_.d"} ) {
     if ( -e $file ) { 
      $Checked{$_} = (stat($file))[9]; $out = check_time ( $i, $file, $time_o )
     }
     last LINE if $out;
    }
   }   
  }
 $out; 
}

#*******************************************************************************
# subroutine read_cpp
#*******************************************************************************
sub read_cpp {

# gets all cpp commands from all include files and adds them to @out
 my $i; $i = $_[0];
 my $file; $file = $_[1];
 my $ext; $ext = " ";
 my $n; $n = 0;
 my @in; undef @in;
 
 if ( -e $file && $Mk{Include_Extension}[$i] ) {
  open IN, $file; while (<IN>) { push ( @in, "$_" ) }; close IN;
  foreach (@in) {
   chomp ( $_ ); s/^ +//;
   if ( $n == 0 ) { $ext = $_; $ext =~ s/.*\.// } 
   if ( $n > 0 ) {
    if ( /^#/ ) { 
     if ( index ( $Mk{Include_Extension}[$i], $ext ) > -1 ) {
      if ( ! /^#\s+include/ && ! /^#include/ ) { push ( @out, "$_\n" ) } 
     } 
    }
    else {read_cpp ( $i, "$Mk{Dependency_Directory}[$i]/$_.d" ) }
   }
   $n++;
  }
 }
 1;
}


#*******************************************************************************
# subroutine mk_dep
#*******************************************************************************
sub mk_dep {

# expects to be in code directory
 my $i; $i = $_[0];
 my $file; $file = $_[1];
 my $file_d; $file_d = $_[2];
 my @dep; my $dep; my @dep_list; my @in; my $dep_list; my $a;
 my $ext; my $base; my @extension_list; my $extension;

 undef @dep; undef @in;
  
 $f = $file; $f =~ s/.*\///g; $f =~ s/\..*//;
 if ( -e $file ) { open IN, $file; while (<IN>) { push ( @in, $_ ) }; close IN }
 foreach (@in) {
  chomp($_); s/USE /use /g; s/INCLUDE /include /g; s/^\s+//g;
  s/^#\s+include /include /g; s/^#include /include /g;
  $dep = $_;

# find include files
  if ( /^include / ) { 
   $dep =~ s/include //g; $dep =~ s/"//g;$dep =~ s/'//g; 
   $dep =~ s/<//g; $dep =~ s/>//g; $dep =~ s/!.*//g;
   if ( $dep ) { push (@dep, $dep) } 
  }

# find module files
  if ( /^use / ) { 
   $dep =~ s/use //g; $dep =~ s/,.*//g; $dep =~ s/!.*//g; $dep =~ s/\s+//g; 
   if ( $dep ) {
    if ( ! %Modules ) { find_mod () }
    else { if ( ! $Modules{$dep} ) { find_mod () } }
    if ( $Modules{$dep} ) {
     $a = $Modules{$dep};  $a =~ s/.*\///g; $a =~ s/\..*//;
#    check the file is not dependent on itself (causing deep recursion)
     if ($a ne $f) {  
      if ( $Mk{Code_Extension}[$i] ) {
       $a = $Modules{$dep};  $a =~ s/.*\///g; $a =~ s/\..*//;
       $a = "$a.$Mk{Code_Extension}[$i]";
       push ( @dep, "$a" );
      }
      if ( $Mk{Object_Extension}[$i] ) {
       $a = $Modules{$dep};  $a =~ s/.*\///g; $a =~ s/\..*//;
       if ($a ne $f) {
        $a = "$a.$Mk{Object_Extension}[$i]";
        $dir_o = $Mk{Object_Directory}[$i]; $dir_o =~ s/.*\///;
        push ( @dep, "$dir_o/$a" );
       }
      }
     }
    }
    else { 
     open OUT, ">>$Mk{Log_File}[$i]";
     print OUT "=>Warning: Can not find module $dep\n";
     close OUT
    }
   }
  }
  
 }
 
 open OUT, ">$file_d"; $a = $file; $a =~ s/.*\///g; print OUT "$a\n"; close OUT;
 $dep_list = $a;
 foreach $dep (@dep) {
 
  $dep =~ s/ +//; $dep =~ s/$Mk{Code_Directory}[$i]//g;
  if ( $dep ) {

   if ( index( $dep_list, $dep ) < 0 ) {
    $ext = $dep;  $ext =~ s/.*\///g; $ext =~ s/.*\.//;
    $base = $file; $base =~ s/.*\///g; $base =~ s/\..*//;
    mk_code ( $dep );
#   make this dependency file depend on each dependencies dependents
    $a = $dep; 
    $a = "$Mk{Dependency_Directory}[$i]/$a.d";
    if ( -e $a ) {
     undef @in; open IN, $a; while (<IN>) { push ( @in, $_ ) }; close IN;
     foreach (@in) {
      chomp($_);
      if (! /^#/ ) {
       $a = $_;
       if ( $a ) {
        if ( index( $dep_list, $a ) < 0 ) {
         mk_code ( $a );
         if ( $a ) {
          if ( -e $a ) {
           open OUT, ">>$file_d"; $a =~ s/.*\///g; print OUT "$a\n"; close OUT;
           $dep_list = "$dep_list $a ";
          }
         }
        }
       }
      }
     }
    }
    
    $a = $dep;
#   make this dependency file depend on each dependency if not already
    if ( $a ) {
     if ( index( $dep_list, $a ) < 0 ) {
      if ( -e $a ) { 
       open OUT, ">>$file_d"; print OUT "$a\n"; close OUT;
       $dep_list = "$dep_list $a ";
      }
     }
    }

   }

  }
 }
 
# add cpp commands from source file
 undef @in;
 $ext = $file;  $ext =~ s/.*\///g; $ext =~ s/.*\.//;
 DIRECTORY:
  for $n ( 0 .. $#{ $Mk{Source_Directory} } ) {
   $file_s = $file; $file_s =~ s/.*\///g;
   $file_s =  "$Mk{Source_Directory}[$n]/$file_s";
   if ( ! -e $file_s ) {
    if ( index( $Mk{Include_Extension}[$n], $ext ) < 0 ) {
     @extension_list = split( / /, $Mk{Source_Extension}[$n] );
     EXTENSION:
      foreach $extension (@extension_list) {
       $file_s = $file; $file_s =~ s/.*\///g; $file_s =~ s/\..*$//;
       $file_s =  "$Mk{Source_Directory}[$n]/$file_s.$extension";
       if ( -e $file_s ) { last EXTENSION }
      }
    }
   }
   if ( $file_s && $Mk{Code_Directory}[$n] && $Mk{Code_Directory}[$i] ) {
    if ( -e $file_s && $Mk{Code_Directory}[$n] eq $Mk{Code_Directory}[$i]) {
     open IN, $file_s; 
     while (<IN>) {
      s/^\s+//; if ( /^#/ ) { push ( @in, $_ ) } 
     }
     close IN; last DIRECTORY;
    }
   }
  }
 open OUT, ">>$file_d"; foreach (@in) { print OUT "$_" }; close OUT;
 undef @dep;
 @dep = split( / /, $dep_list );
 foreach $d (@dep) {
  if ( "$d" ne "" ) {
   for $n ( 0 .. $#{ $Mk{Source_Directory} } ) {
    $d = "$Mk{Dependency_Directory}[$n]/$d.d";
    if ( -e $d && $Mk{Code_Directory}[$n] eq $Mk{Code_Directory}[$i]) {
     undef @in; $n = 0;
     open IN, "$d"; 
     while (<IN>) {
      s/^\s+//; if ( /^#/ ) { 
       push ( @in, $_ ); s/\#//; s/^\s+//; if ( /^define / ) { $n++ }
      }
     }
     close IN;
     if ( $n > 0 ) {
      open OUT, ">>$file_d"; foreach (@in) { print OUT "$_" }; close OUT;
     }
    }
   }
  }
 }
# open OUT, ">>$file_d"; foreach (@in) { print OUT "$_" }; close OUT;

 1;
}

#*******************************************************************************
# subroutine find_mod
#*******************************************************************************
sub find_mod {

# finds files that contain modules and puts the file name in %Modules
# global: %Modules, %Mk, $idx
 my $i; my @extension_list; my $extension; my @file_list; my $file; 

# search all source directories for modules
 for $i ( 0 .. $#{ $Mk{Source_Directory} } ) {
# check if this source directory shares a code directory with the target
  if ( $Mk{Code_Directory}[$i] eq $Mk{Code_Directory}[$idx] ) {
#  get all source extensions for this source directory
   @extension_list = split( / /, $Mk{Source_Extension}[$i] );
   foreach $extension (@extension_list) {
    @file_list = <$Mk{Source_Directory}[$i]/*.$extension>;
#   search all source files that have a Source_Extension
    foreach $file (@file_list) {
     open IN, $file;
#    read file
     while (<IN>) {
      chomp ($_);
      s/ //g; s/MODULEPROCEDURE//g; s/moduleprocedure//g; s/MODULE/module/g;
#     set %Modules if a module is found
      if ( /^module/ ) { s/module//g; $file =~ s/.*\///; $Modules{$_} = $file }
     }
     close IN;
    }
   }
  }
 }
 1;
}

#*******************************************************************************
# subroutine set_var
#*******************************************************************************
sub set_var {

# adds path $dir if path is not already set
 my $var; my $dir;

 $var = $_[0];
 $dir = $_[1];
 if ( $var ) {
  $var =~ s/ //g;
  if ( $dir ) {
   $dir =~ s/ //g;
   if ( index ( $var, "./") == 0 ) { $var =~ s/\.\///; $var = "$dir/$var" }
   if ( index ( $var, "/") != 0 ) { $var = "$dir/$var" }
  }
  $var =~ s/\/\//\//g; $var =~ s/ //g;
 }
 $var;
}

#*******************************************************************************
# subroutine read_mk
#*******************************************************************************
sub read_mk {

# reads a "mk.in" file to set mk environment variables (reads into %MkH)
# global: %MkH, $Definitions_Up, $Definitions_Down, $Equality_Character

 local $file; local $num; local $var; local $idx; local $val;
 $file = $_[0];

# set some variables if not already set
 unless ( $Definitions_Up ) { $Definitions_Up = "true" }
 unless ( $Definitions_Down ) { $Definitions_Down = "false" } 
 unless ( $Equality_Character ) { $Equality_Character = "=" }

 $Equality_Character =~ s/^\s+//; $Equality_Character =~ s/ .*//;

# set $num to the number of source directories before this file is read.
 $num = $#{ $MkH{Source_Directory} };

# read the file if it exists
 if ( -e $file ) {

  open ( MkH_IN, $file );
  
  LINE: 
   while ( <MkH_IN> ) {

#   remove extra blanks and any comments (started with #)
    chomp($_); s/^\s+//; s/\s+$//; s/#.*//;
    
#   check for changed $Equality_Character 
    if ( /Equality_Character/ ) {
     s/.*Equality_Character//; s/^\s+//; s/ .*//; 
     $Equality_Character = $_; $_ = "";
    }
    $val = $_;
#   split on first $Equality_Character
    if ( /$Equality_Character/ ) { 
     ( $var, $val ) = split /$Equality_Character/, $_, 2;
#    replace other brackets with parethesis
     $var =~s/\[/\(/g; $var =~s/\]/\)/g; $var =~s/\{/\(/g; $var =~s/\}/\)/g;
     $var =~ s/ //g;
#    look for an array index
     $idx = $var; $idx =~ s/\).*//; $idx =~ s/.*\(//;
     if ( $idx eq $var ) { $idx = 0 }
      $var =~ s/\(.*//;
#    if variable is already defined skip to the next variable 
     if ( $MkH{$var}[$idx] ) { undef $var; next LINE }
#    if sources are already defined by another file skip to the next variable 
     if ( $var eq "Source_Directory" and $num >= 0 ) { undef $var; next LINE }
    }
    if ( $var && $val ) {
#    add to the value of this variable
     $val =~ s/^\s+//; $val =~ s/\s+$//; $val =~ s/'//g; $val =~ s/"//g;
     if ( $MkH{$var}[$idx] ) { $MkH{$var}[$idx] = "$MkH{$var}[$idx] $val" }
     else { $MkH{$var}[$idx] = $val }
    }

   }
  
  close MkH_IN;
# if a variable is not defined, get definition from previous source definition
  foreach $var ( keys %MkH ) {
   for ( $idx = 1; $idx <= $#{ $MkH{Source_Directory} }; $idx++ ) {
    unless ( $MkH{Definitions_Up}[$idx] ) { 
     $MkH{Definitions_Up}[$idx] = $Definitions_Up;
    }
    unless ( $MkH{Definitions_Down}[$idx] ) {
     $MkH{Definitions_Down}[$idx] = $Definitions_Down;
    }
    if ( $MkH{Definitions_Up}[$idx] eq "true" ) {
     if ( $var ne "Source_Directory" ) {
      unless ( $MkH{$var}[$idx] ) {$MkH{$var}[$idx] = $MkH{$var}[$idx-1] }
     }
    }
   }
  }
  
# if a variable is not defined, get definition from next source definition
  foreach $var ( keys %MkH ) {
   for ( $idx = $#{ $MkH{Source_Directory} }-1; $idx >= 0; $idx-- ) {
    unless ( $MkH{Definitions_Up}[$idx] ) { 
     $MkH{Definitions_Up}[$idx] = $Definitions_Up;
    }
    unless ( $MkH{Definitions_Down}[$idx] ) {
     $MkH{Definitions_Down}[$idx] = $Definitions_Down;
    }
    if ( $MkH{Definitions_Down}[$idx] eq "true" ) {
     if ( $var ne "Source_Directory" ) {
      unless ( $MkH{$var}[$idx] ) {$MkH{$var}[$idx] = $MkH{$var}[$idx+1] }
     }
    }
   }
  }
  
# eliminate definitions without associated source directories
  $num = 0;
  for $idx ( 0 .. $#{ $MkH{Source_Directory} } ) {
   if ( $MkH{Source_Directory}[$idx] ) {
    foreach $var ( keys %MkH ) {
     $MkH{$var}[$num] = $MkH{$var}[$idx];
    }
    $num++;
   }
  }  
  $#{ $MkH{Source_Directory} } = $num - 1;

 }   
 1;
}


#*******************************************************************************
# subroutine find_path
#*******************************************************************************
sub find_path {

# finds a file given a search path and adds path to the file
 my $file; my $path_search; my $path_mk; my $path; my @path_list; my $path_file;
 $file = $_[0];
 $path_search = $_[1]; unless ( $path_search ) { $path_search = " " }
 $path_mk = $_[2]; unless ( $path_mk ) { $path_mk = " " }
 
 if ( $file ) {
  $path_file = $file;
# create the full search path (file + environment)
  $path = "$path_search $path_mk"; $path =~ s/:/ /g;
  @path_list = split( / /, $path );
  PATH:
   foreach $path (@path_list) {
    $path_file = "$path/$file"; $path_file =~ s/ //g;
    if ( -e "$path_file") { last PATH }
   }
 }
 $path_file;
}


#*******************************************************************************
# subroutine mk_unique
#*******************************************************************************
sub mk_unique {

# removes adjacent similar items from an array
 my @items; my $item;
 @items=" ";
 
 foreach $item (@_) { 
  if ( $item ne $items[$#items] ) { push ( @items, $item ) } 
 }
 shift(@items);
 @items;
}


#*******************************************************************************
# subroutine write_out
#*******************************************************************************
sub write_out {

# writes an array of words to a certain column width to filehandle OUT
 my $column; my $tab; my $width;
 $column = 0;
 $tab = "   ";
 $width = 80;
 
 print OUT $tab; $column += length($tab);
 foreach (@_) {
  s/^ +//; s/$ +//;
  $column += length(" $_");
  if ($column > $width) {print OUT "\n$tab"; $column = length("$tab $_")}
  print OUT " $_";
 }
 print OUT "\n";
 1;
}

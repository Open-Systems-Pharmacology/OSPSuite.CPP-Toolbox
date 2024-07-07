@echo off

nuget pack OSPSuite.CPP-Toolbox_Windows.nuspec -version %1
IF %ERRORLEVEL% NEQ 0 EXIT /B %ERRORLEVEL%

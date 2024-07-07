@echo off
setlocal

if not exist Build mkdir Build
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"
IF %ERRORLEVEL% NEQ 0 EXIT /B %ERRORLEVEL%

rem for %%T in (Debug RelWithDebInfo) do (
for %%T in (Debug Release) do (
    echo Compiling for build type = %%T   
    cmake -BBuild/%%T/x64/ -Hsrc/OSPSuite.CPP-Toolbox/ -DCMAKE_BUILD_TYPE=%%T -DBUILD_SHARED_LIBS=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON 
    IF %ERRORLEVEL% NEQ 0 EXIT /B %ERRORLEVEL%
    msbuild Build/%%T/x64/ALL_BUILD.vcxproj -p:Configuration=%%T
    IF %ERRORLEVEL% NEQ 0 EXIT /B %ERRORLEVEL%
    if exist Dist\Windows\%%T\x64 rmdir /S /Q Dist\Windows\%%T\x64
    IF %ERRORLEVEL% NEQ 0 EXIT /B %ERRORLEVEL%
    mkdir Dist\Windows\%%T\x64
    IF %ERRORLEVEL% NEQ 0 EXIT /B %ERRORLEVEL%
)

for /R Build\Debug\x64\Debug %%F in (*.dll *.lib *.pdb) do copy "%%F" Dist\Windows\Debug\x64\
IF %ERRORLEVEL% NEQ 0 EXIT /B %ERRORLEVEL%
for /R Build\Release\x64\Release %%F in (*.dll *.lib *.pdb) do copy "%%F" Dist\Windows\Release\x64\
IF %ERRORLEVEL% NEQ 0 EXIT /B %ERRORLEVEL%

endlocal
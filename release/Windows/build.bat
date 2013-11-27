msbuild.exe idock.vcxproj /t:Rebuild /p:Configuration=Release;Platform=x64
"%CUDA_PATH%\bin\nvcc" -arch=compute_11 -fatbin -ccbin "%VS110COMNTOOLS%\..\..\VC\bin\amd64" -o ..\..\bin\idock.fatbin ..\..\src\idock.cu

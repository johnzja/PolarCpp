mkdir('Release');
disp('Compiling...');
mex -R2018a -DON_LINUX -silent Qary_SC_Decoder/dllmain.cpp Qary_SC_Decoder/pch.cpp PolarCpp/GF.cpp PolarCpp/Qary_dist.cpp PolarCpp/SC.cpp PolarCpp/SCFrame.cpp PolarCpp/SCList.cpp ...
    -output Release/Qary_SC_Decoder
disp('MEX successfully generated at ./Release/');
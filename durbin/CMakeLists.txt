xilinx_synthesis_target("Durbin" "" "Durbin" "")
xilinx_synthesis_target("DurbinMod" "" "DurbinMod" "")
xilinx_testbench("TestMod.cpp DurbinMod.cpp" "DurbinMod")
xilinx_testbench("Test.cpp Durbin.cpp" "Durbin")

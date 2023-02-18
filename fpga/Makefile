VPP := $(XILINX_VITIS)/bin/v++
EMCONFIGUTIL := $(XILINX_VITIS)/bin/emconfigutil
TARGET:= sw_emu
PLATFORM := xilinx_u250_gen3x16_xdma_2_1_202010_1

# kernel targets
KRNL_SCALE_XO := krnl_scale.$(TARGET).xo
#KRNL_VMULT_XO := krnl_vmult.$(TARGET).xo
XCLBIN := kernel.$(TARGET).xclbin

# host target
HOST_SRC := ./host/main.cpp ./host/utils.cpp 
HOST_EXE := host.exe

# config files target
EMCONFIG_FILE := emconfig.json

VPP_OPTS := -s -t $(TARGET) --platform $(PLATFORM)
VPP_OPTS += --temp_dir ./all_builds/$(TARGET)
VPP_OPTS += --report_dir ./all_logs/$(TARGET)
VPP_OPTS += --log_dir ./all_logs/$(TARGET)

CFLAGS := -g -std=c++1y -I$(XILINX_XRT)/include 
#LFLAGS := -L$(XILINX_XRT)/lib -lOpenCL -lxrt_coreutil -lrt -luuid
LFLAGS := -L$(XILINX_XRT)/lib -lxilinxopencl -pthread -lrt -lxrt_coreutil
NUMDEVICES := 1

# run time args
CMD_ARGS:= -f kernel.${TARGET}.xclbin

# primary build targets
.PHONY: xclbin host all

xclbin:  $(XCLBIN)
host: $(HOST_EXE) 

all: xclbin host 

clean:
	-$(RM) $(EMCONFIG_FILE) $(HOST_EXE) $(XCLBIN) 

cleanall: clean	
	rm -rf all_* *.log *.xclbin *.compile_summary *.link_summary *.info *.xo

# kernel rules
$(KRNL_SCALE_XO): ./kernels/scale/scale.cpp
	$(RM) $@
	$(VPP) $(VPP_OPTS) -c -k scale_kernel -I'$(<D)' -o $@ $+

#$(KRNL_VMULT_XO): ./vmult_kernel/krnl_vmult.cpp
#	$(RM) $@
#	$(VPP) $(VPP_OPTS) -c -k krnl_vmult -I'$(<D)' -o $@ $+

$(XCLBIN): $(KRNL_SCALE_XO)
#$(XCLBIN): $(KRNL_VADD_XO) $(KRNL_VMULT_XO) 
	$(VPP) $(VPP_OPTS) -l -o $@ $(KRNL_XO) $+ 

# host rules
$(HOST_EXE): $(HOST_SRC) 
	g++ $(CFLAGS) -o $@ $+ $(LFLAGS)
	@echo 'Compiled Host Executable: $(HOST_EXE)'

$(EMCONFIG_FILE):
	$(EMCONFIGUTIL) --nd $(NUMDEVICES) --od . --platform $(PLATFORM)

run: $(XCLBIN) $(HOST_EXE) $(EMCONFIG_FILE)
	XCL_EMULATION_MODE=${TARGET} ./$(HOST_EXE) $(CMD_ARGS)

.PHONY: help

help::
        $(ECHO) "Makefile Usage:"
	$(ECHO) "  make all TARGET=<sw_emu/hw_emu/hw> PLATFORM=<FPGA platform> LAB=<run1/run2>"
	$(ECHO) "      Command to generate the design for specified Target and Device for LAB step 1 or 2."
	$(ECHO) ""
	$(ECHO) "  make xclbin TARGET=<sw_emu/hw_emu/hw> PLATFORM=<FPGA platform> LAB=<run1/run2>"
	$(ECHO) "      Command compile just the kernel of the design for specified Target and Device for LAB step 1 or 2."
	$(ECHO) ""
	$(ECHO) "  make host TARGET=<sw_emu/hw_emu/hw> PLATFORM=<FPGA platform>"
	$(ECHO) "      Command to generate the host application for specified Target and Device."
	$(ECHO) ""
	$(ECHO) "  make run TARGET=<sw_emu/hw_emu/hw> PLATFORM=<FPGA platform> LAB=<run1/run2>"
	$(ECHO) "      Command compile just the kernel of the design for specified Target and Device for LAB step 1 or 2."
	$(ECHO) ""
	$(ECHO) "  make clean "
	$(ECHO) "      Command to remove the generated non-hardware files."
	$(ECHO) ""
	$(ECHO) "  make cleanall"
	$(ECHO) "      Command to remove all the generated files."
	$(ECHO) ""

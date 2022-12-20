SM=52
ARCH = -arch sm_$(SM)

source = RKF.cu Host.cu perturbers.cu read.cu output.cu interpolate.cu Scan.cu
headers = define.h Host.h force.h integrator.h interpolateGPU.h

RKF: $(source) $(headers)
	nvcc $(ARCH) --compiler-options -Wall -o RKF $(source) $(BUILD_DATA)


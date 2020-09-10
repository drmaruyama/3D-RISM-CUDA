CUDA_DIR=/usr/local/cuda

CC 	= g++

CFLAGS 	= -O3 -fomit-frame-pointer -fopenmp

NVCC	= nvcc

NVFLAGS = -O3 -m64 -gencode=arch=compute_70,code=sm_70 -I$(CUDA_DIR)/include 

LDFLAGS = -L$(CUDA_DIR)/lib64 
LIBS = -lcufft -lcuda -lcudart 

OBJS   	= main.o\
	  add_tuv.o\
	  alloc2D.o\
	  alloc3D.o\
	  anderson.o\
	  anderson_cuda.o\
	  cal_Coulomb.o\
	  cal_euv.o\
	  cal_exchem.o\
	  cal_exdft.o\
	  cal_grad.o\
	  cal_LJ.o\
	  cal_pmv.o\
	  cal_potential.o\
	  cal_rms.o\
	  calculate.o\
	  calculate_cuda.o\
	  calloc2D.o\
	  cell.o\
	  dealloc2D.o\
	  dealloc3D.o\
	  excp_fd.o\
	  fft3d.o\
	  initialize.o\
	  initialize_g.o\
	  initialize_tuv.o\
	  iterate.o\
	  output.o\
	  output_cuv.o\
	  output_euv.o\
	  output_grad.o\
	  output_guv.o\
	  output_huv.o\
	  output_xmu.o\
	  read_input.o\
	  read_tuv.o\
	  set_cuda.o\
	  set_fname.o\
	  set_solvent.o\
	  solute.o\
	  solvent_read.o\
	  solvent_spline.o\
	  solvent_spline2.o\
	  spline.o\
	  splint.o\
	  write_tuv.o\
	  
PROGRAM = 3drism-cuda

all:	$(PROGRAM)

$(PROGRAM) : $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(NVFLAGS) -c $< -o $@

clean:;	@rm -f $(OBJS) core

$(OBJS) : 

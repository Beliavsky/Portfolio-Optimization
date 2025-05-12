executables = xportfolio_opt_sim_gfort.exe
FC     = gfortran
FFLAGS = -O0 -Wall -Werror=unused-parameter -Werror=unused-variable -Werror=unused-function -Wno-maybe-uninitialized -Wno-surprising -fbounds-check -static -g -fmodule-private
obj    = kind.o portfolio_opt.o xportfolio_opt_sim.o

all: $(executables)

# Compile .f90 to .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

xportfolio_opt_sim_gfort.exe: kind.o portfolio_opt.o xportfolio_opt_sim.o
	$(FC) -o xportfolio_opt_sim_gfort.exe kind.o portfolio_opt.o xportfolio_opt_sim.o $(FFLAGS)

run: $(executables)
	./xportfolio_opt_sim_gfort.exe

clean:
	rm -f $(executables) $(obj)


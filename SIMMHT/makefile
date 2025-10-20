# =========================
# Makefile (gfortran)
# =========================

# Compilador
FC := gfortran

# Pastas de build
OBJDIR := obj
MODDIR := mod

# Fontes (na ordem correta p/ módulos)
SRCS := variables.f90 functions.f90 input.f90 simmht.f90
OBJS := $(SRCS:%.f90=$(OBJDIR)/%.o)

# Flags
# -ffree-line-length-none  evita truncamento de linhas longas (seus arquivos usam linhas > 132)
# -fopenmp                 OpenMP
# -march=native -O3        otimizações (release)
# -g -fcheck=all ...       debug amigável
FFLAGS_BASE    := -fopenmp -J$(MODDIR) -ffree-line-length-none
FFLAGS_RELEASE := -O3 -march=native -funroll-loops $(FFLAGS_BASE)
FFLAGS_DEBUG   := -O0 -g -fcheck=all -fbacktrace -Wall -Wextra $(FFLAGS_BASE)

# Alvos
.PHONY: all release debug clean

all: release

release: FFLAGS=$(FFLAGS_RELEASE)
release: simmht.ex

debug: FFLAGS=$(FFLAGS_DEBUG)
debug: simmht_debug.ex

# Linkagem
simmht.ex: $(OBJS) | $(MODDIR)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

simmht_debug.ex: $(OBJS) | $(MODDIR)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

# Regra de compilação
$(OBJDIR)/%.o: %.f90 | $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Pastas
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(MODDIR):
	mkdir -p $(MODDIR)

# Limpeza
clean:
	rm -rf $(OBJDIR) $(MODDIR) simmht.ex simmht_debug.ex

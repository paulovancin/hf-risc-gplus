CC = gcc
CFLAGS = -O2 -Wall

all:
	$(CC) $(CFLAGS) hf_riscv_sim.c -o hf_riscv_sim
	$(CC) $(CFLAGS) hf_riscv_sim_full.c -o hf_riscv_sim_full

clean:
	rm -rf *.o *~ hf_riscv_sim hf_riscv_sim_full

baselineAndEfficiency:baselineAndEfficiency.c
	gcc -lm -lpthread -o baselineAndEfficiency baselineAndEfficiency.c
clean:
	rm baselineAndEfficiency
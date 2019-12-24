# Genetic algorithm
This python script uses a genetic algorithm to find the maximum value of a **quadratic equation**.

20 chromosomes are used through 50 generations. The value for each chromosome is **encoded in a 22 bit number**. The higher values have a better chance for **crossover**. Also, in this approach, **mutation** is flipping a random bit from a chromosome.

> **Note:**  The quadratic equation, number of chromosomes/generations and crossover/mutation probabilities can be updated from the input file.

### The script displays a plot with the evolution of the maximum value through the generations
![Example](https://github.com/andrei828/GeneticAlgorithm/blob/master/plot.png)
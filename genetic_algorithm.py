import copy
import matplotlib.pyplot as plt
from random import random, randint

DISPLAY_GENERATION = 0
generation = 0

def generate_binary_str(len):
		binary = []
		for i in range(len):
			if random() <= 0.5:
				binary.append('0')
			else:
				binary.append('1')
		return ''.join(binary)

class Chromozome:

	def __init__(self, length=None):
		self.length = length
		self._binary = generate_binary_str(self.length)

	@property
	def binary(self):
		return self._binary
	
	@binary.setter
	def binary(self, value):
		self._binary = value

	@property
	def selection_prob(self):
		return self._selection_prob

	@selection_prob.setter
	def selection_prob(self, value):
		self._selection_prob = value
	
	def decode_binary(self, func_domain, precision):
		decimal = 0
		for gene_index in range(self.length):
			if self.binary[self.length - gene_index - 1] == '1':
				decimal += 2**gene_index
		
		return round (
			decimal * (func_domain[1] - func_domain[0]) / 
			(2**self.length - 1) + func_domain[0], precision )


def print_population(population, func_domain, func_coef, precision):
	for i, chromozome in enumerate(population, start=0):
		print( '{}:\t'.format(i + 1), 
			   chromozome.binary, '\tx=', 
			   chromozome.decode_binary(func_domain, precision), '\tf=', 
			   get_f(chromozome.decode_binary(func_domain, precision), 
			   		 func_coef) )

def print_status(operation, population, func_domain, func_coef, precision):
	if generation == DISPLAY_GENERATION:
		print("After {}".format(operation))
		print_population(population, func_domain, func_coef, precision)

def binary_search(arr, start, finish, element):

	if (start <= finish):
		mid = start + (finish - start) // 2

		if element >= arr[mid] and element < arr[mid + 1]:
			return mid + 1
		elif element > arr[mid]:
			return binary_search(arr, mid + 1, finish, element)
		else:
			return binary_search(arr, start, mid - 1, element)

	return -1

def get_chromozome_length(func_domain, precision):
	val = (func_domain[1] - func_domain[0]) * 10**precision

	chrm_length = 1
	while 2**chrm_length < val:
		chrm_length += 1
	return chrm_length

def get_f(x, func_coef):
	return (1.0 * func_coef[0] * x**2 + func_coef[1] * x + func_coef[2])

def get_max_value(population, func_domain, func_coef, precision):
	return max([get_f(chromozome.decode_binary(func_domain, precision), 
				func_coef) for chromozome in population] )

def get_average_values(population, func_domain, func_coef, precision):
	sum = 0
	for chromozome in population:
		sum += get_f(chromozome.decode_binary(func_domain, precision), func_coef)
	return sum / len(population)

def set_selection_prob(population, func_domain, precision, func_coef):
	
	if generation == DISPLAY_GENERATION:
		print("Selection probability")

	for i, chromozome in enumerate(population, start=0):

		x = chromozome.decode_binary(func_domain, precision)
		f_x = get_f(x, func_coef)
		
		sum = 0.0
		for chr in population:
			x_chr = chr.decode_binary(func_domain, precision)
			f_x_chr = get_f(x_chr, func_coef)
			sum += f_x_chr

		chromozome.selection_prob = f_x / sum
		if generation == DISPLAY_GENERATION:
			print('chromozome ', i + 1, ' probability ', chromozome.selection_prob)

	selection_interval = [0, population[0].selection_prob]
	for chromozome in population[1:]:
		selection_interval.append(selection_interval[-1] + chromozome.selection_prob)

	if generation == DISPLAY_GENERATION:
		print("Selection interval")
		for interval in selection_interval:
			print(interval, end=' ')
		print()

	return selection_interval

def crossover(chromozome_1, chromozome_2):
	crossover_barrier = randint(0, chromozome_1.length - 1)
	binary_1 = list(chromozome_1.binary)
	binary_2 = list(chromozome_2.binary)
	
	if generation == DISPLAY_GENERATION:
		print( chromozome_1.binary, chromozome_2.binary, 
			   'barrier ', crossover_barrier)

	for i in range(crossover_barrier):
		tmp = binary_1[i]
		binary_1[i] = binary_2[i]
		binary_2[i] = tmp

	chromozome_1.binary = ''.join(binary_1)
	chromozome_2.binary = ''.join(binary_2)

	if generation == DISPLAY_GENERATION:
		print("Result\t", chromozome_1.binary, chromozome_2.binary)

def mutation(chromozome):
	binary = list(chromozome.binary)
	num_of_bit_flips = 1#randint(0, chromozome.length - 1)

	while num_of_bit_flips:
		index = randint(0, chromozome.length - 1)
		binary[index] = '1' if binary[index] == '0' else '0'
		num_of_bit_flips -= 1

	chromozome.binary =  ''.join(binary)

def crossover_selection(population, recombine_prob):
	
	if generation == DISPLAY_GENERATION:
		print("Crossover probability: ", recombine_prob)

	selected_chromozomes = []
	selected_chromozomes_index = []
	for i, chromozome in enumerate(population, start=0):
		probability = random()

		if generation == DISPLAY_GENERATION:
			print('{}:\t'.format(i + 1), chromozome.binary, '\tu=', probability, end='')

		if probability < recombine_prob:
			selected_chromozomes.append(chromozome)
			selected_chromozomes_index.append(i + 1)
			if generation == DISPLAY_GENERATION:
				print(' <', recombine_prob, end='')
		
		if generation == DISPLAY_GENERATION:
			print()
	
	if len(selected_chromozomes) % 2:
		selected_chromozomes = selected_chromozomes[:-1]
	
	for i in range(0, len(selected_chromozomes) - 1, 2):
		if generation == DISPLAY_GENERATION:
			print("Crossover between chromozome", selected_chromozomes_index[i],
				  "and chromozome ", selected_chromozomes_index[i + 1])

		crossover(selected_chromozomes[i], selected_chromozomes[i + 1])

def mutation_selection(population, mutation_prob):

	no_mutation = True
	if generation == DISPLAY_GENERATION:
		print("Mutation probability: ", mutation_prob)
		print("Following chromozomes suffered mutation")

	for i, chromozome in enumerate(population, start=0):
		probability = random()

		if probability < mutation_prob:
			if generation == DISPLAY_GENERATION:
				print(i + 1)
			no_mutation = False
			mutation(chromozome)

	if generation == DISPLAY_GENERATION and no_mutation:
		print("Not a single one")


def rearrange_chromozomes_in_population(selection_interval, population, func_domain, precision, func_coef):
	new_population = []
	
	for index in range(len(population)):
		probability = random()

		selected_index = binary_search( selection_interval, 
										0, len(population) - 1,
										probability )

		if generation == DISPLAY_GENERATION:
			print("u=", probability, "\tselect chromozome ", selected_index)
		new_population.append(copy.deepcopy(population[selected_index - 1]))

	# del population
	return new_population

# read input and create variables
population = []
with open('input_genetic_algorithm.txt', 'r') as document:
	func_coef = [int(coef) for coef in next(document).split()]
	func_domain = [int(margin) for margin in next(document).split()]
	population_size = int(next(document).split()[0])
	precision = int(next(document).split()[0])
	recombine_prob = float(next(document).split()[0])
	mutation_prob = float(next(document).split()[0])
	max_num_generations = int(next(document).split()[0])

	chromozome_length = get_chromozome_length(func_domain, precision)
	for chromozome in range(population_size):
		full_line = next(document).split()
		c = Chromozome(chromozome_length)
		# c.binary = full_line[1]
		population.append(c)

max_values = []
average_values = []
print_population(population, func_domain, func_coef, precision)

while generation < max_num_generations:

	# store important data for each generation
	max_values.append(get_max_value(population, func_domain, func_coef, precision))
	average_values.append(get_average_values(population, func_domain, func_coef, precision))

	selection_interval = set_selection_prob(population, func_domain, precision, func_coef)
	population = rearrange_chromozomes_in_population(selection_interval, population, func_domain, precision, func_coef)
	print_status("selection", population, func_domain, func_coef, precision)

	crossover_selection(population, recombine_prob)
	print_status("crossover", population, func_domain, func_coef, precision)

	mutation_selection(population, mutation_prob)
	print_status("mutation", population, func_domain, func_coef, precision)

	generation += 1

max_values.append(get_max_value(population, func_domain, func_coef, precision))
average_values.append(get_max_value(population, func_domain, func_coef, precision))

print("Max value evolution")
for max_value in max_values:
	print(max_value)

plt.plot(range(max_num_generations + 1), average_values, 'bo')
plt.margins(0.3)
plt.ylabel('average value')
plt.xlabel('generation')
plt.show()








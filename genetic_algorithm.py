import random
import subprocess
import math

POPULATION_SIZE = 30         # Initial population size. This should stay relatively constant.
ITERATIONS = 30              # Amount of generations to go through
PROPORTION_OFFSPRING = 0.8   # 80:20 for offspring:old generation
MUTATION_PROBABILITY = 1000   # 1/x chance of mutation



class Organism(object):
    def __init__(self, bitstring):
        self.bitstring = bitstring
        self.JSONObject = dictToJSON(self.bitstringToDict(0))
        self.fitness = self.getFitness()

    def getFitness(self):
        """Return the fitness of this organism.
        This function should be catered to whatever task you're solving.
        """
        return fitness_distanceMovedPerTime(self)

    def breed(self, partner):
        """Breed a child with 'partner'.

        Uses the cross-over method, as well as applying all necessary
        random mutations to the child.P
        """
        child = self.crossover(partner)
        child.mutate()

        return child
    
    def crossover(self, partner):
        """Apply the cross-over reproduction method and return one of the children.

        Randomly choose a split-point, then create two children based on this point.
        The first child is the male[:splitpoint] + female[splitpoint:]
        The second is female[:splitpoint] + male[splitpoint:]
        
        The function then randomly picks one of these children and returns them.
        Currently, this only does a single split point. It may be prudent to
        change this method to allow multiple splits.
        """
        splitpoint = random.randint(1, min([len(self.bitstring), len(partner.bitstring)]))

        newBitstring = random.choice([
                self.bitstring[:splitpoint] + partner.bitstring[splitpoint:],
                partner.bitstring[:splitpoint] + self.bitstring[splitpoint:]])
                
        newChild = Organism(newBitstring)
        return newChild
        
    def mutate(self):
        """Occasionally flip certain bits in the bitstring.
        The regularity of this is denoted by MUTATION_PROBABILITY.
        """
        for i in range(len(self.bitstring)):
            if blueMoon(): # Once in a blue moon!
                self.bitstring = self.bitstring[:i] + flip(self.bitstring[i]) + self.bitstring[i+1:]

    def bitstringToDict(self, offset):
        """Convert the organism's bitstring to a JSON object."""
        chunkSize = 24

        strNumConnections = self.bitstring[offset:3+offset]
        strAngle = self.bitstring[3+offset:8+offset]
        strLeng = self.bitstring[8+offset:12+offset]
        strFreq = self.bitstring[12+offset:16+offset]
        strAmpl = self.bitstring[16+offset:20+offset]
        strPhase = self.bitstring[20+offset:24+offset]

        numConnections = int(strNumConnections, 2)
        
        limb = {}
        limb["angle"] = (int(strAngle, 2) / 31.0) * math.pi*2.0
        limb["length"] = int(strLeng, 2) * 5.0
        limb["frequency"] = int(strFreq, 2) / 4.0
        limb["amplitude"] = (int(strAmpl, 2) / 15.0) * math.pi
        limb["phase"] = (int(strPhase, 2) / 15.0) * math.pi*2.0
        limb["treeSize"] = 1
        limb["connections"] = []
        
        offset += chunkSize
        
        for i in range(numConnections):
            if offset >= len(self.bitstring):
                break
                
            child = self.bitstringToDict(offset)
            limb["connections"].append( child )
            limb["treeSize"] += child["treeSize"]
            offset += limb["treeSize"]*chunkSize
        
        return limb

    def __cmp__(self, other):
        """A comparison function for the sorting of the population.
        Compare organisms based on their fitness.
        """
        return -cmp(self.fitness, other.fitness)

    def __add__(self, partner):
        """Breed two organisms using the '+' operator.
        
        This is in no way necessary -- I just wanted to take advantage of Python's 
        operator-overloading!
        """
        return self.breed(partner)
           


def main():
    print 'Program has started!'
    population = makePopulation()
  #  print 'Population created successfully.'
  #  print [x.fitness for x in population]
    for i in xrange(ITERATIONS):
        population = nextGeneration(population)
        print 'Average population fitness at generation %d (%d%%): %f' % (i+1, 100 * i/ITERATIONS, averagePopulationFitness(population))


    raw_input('Optimisation finished. Press enter to view your creations!')

    for h in sorted(population):
        print 'aoeu'
        viewSimulation(h)

    raw_input('Program complete. Press enter to exit.')
    

def nextGeneration(population):
    """Make the next generation of organisms, given the current generation.

    First, sort the population in terms of descending fitness (the Organism class 
    already has a __cmp__() method which allows us to do this easily).

    Next, reserve the best (POPULATION_SIZE * (1 - PROPORTION_OFFSPRING)) members
    of our population for the next generation. This will guarantee that our 
    algorithm continues to climb hills.

    Breed our new generation of (POPULATION_SIZE * PROPORTION_OFFSPRING) organisms.
    The parents are randomly chosen from the set of old-gen organisms with a
    fitness greater than the average fitness of that population.

    TODO: Add roulette-wheel selection (fitness proportionate selection).
    """
    population.sort()
    newPopulation = population[:int(POPULATION_SIZE * (1 - PROPORTION_OFFSPRING))]

    culledPopulation = applyNaturalSelection(population)
  #  print culledPopulation
    i = 0

    while i < POPULATION_SIZE * PROPORTION_OFFSPRING:
        father = random.choice(culledPopulation)
        mother = random.choice(culledPopulation)

        child = father + mother  # W00T for operator-overloading!!!
        while child.fitness is None:
            father = random.choice(culledPopulation)
            mother = random.choice(culledPopulation)
            child = father + mother

        newPopulation.append(child)
        i += 1
        
    
    return newPopulation
        

def makePopulation():
    """Create a population of POPULATION_SIZE random organisms.

    Note that this function does not actually create the organisms, it calls
    on the generateOrganismBitstring() function for this.
    """
    population = []

    for i in range(POPULATION_SIZE):
        population.append(generateOrganism())
        print 'Population %d%% generated' % ((i*100)/POPULATION_SIZE)
    return population
    

def averagePopulationFitness(population):
    """Get the average fitness of a population."""
    return sum([org.fitness for org in population]) / float(len(population))


def medianPopulationFitness(population):
    """Same as the above but for median, not average."""
    return population[len(population)/2].fitness


def applyNaturalSelection(population, fitnessDistribution=averagePopulationFitness):
    """Kill off all organisms with a fitness level below the average fitness
    of the population.

    The kwarg 'fitnessDistribution' is the function to find the average/median/etc
    """
    averageFitness = fitnessDistribution(population)
    culledPopulation = filter(lambda org: org.fitness >= averageFitness, population)

    return culledPopulation


def generateOrganism():
    """Generate an organism.
    This function should be catered to whatever task you're solving.
    """
    return Organism(makeLimb(MAX_LIMBS))


def blueMoon():
    """Has a 1/MUTATION_PROBABILITY chance of returning true.
    Useful for deciding when to mutate.
    """
    return random.randint(0, MUTATION_PROBABILITY) == 0


def flip(bit):
    """Flip a single bit in a bitstring. '1' -> '0' and vice versa."""
    return '1' if bit == '0' else '0'


#############################
## Humperdink-specific code
#############################

MAX_LIMBS     = 16
MAX_FREQUENCY = 12
MAX_ANGLE     = 31
MAX_LENGTH    = 15
MAX_AMPLITUDE = 15
MAX_PHASE     = 15
MAX_CHILDREN  = 8

def decimalToBinary(n):
    if n == 0: return ''
    return decimalToBinary(n >> 1) + str(n & 1)


def paddedBinary(n, size):
    """Convert a decimal number to binary, and pad the number at the start with 0's"""
    pure = decimalToBinary(n)
    return ('0' * (size - len(pure))) + pure


def dictToJSON(d):
    return str(d).replace("'", '"')


def makeLimb(limbsPermitted, conv=paddedBinary):
    """Make a limb and any child limbs that go with it.

    'limbsPermitted' denotes the amount of limbs/child limbs the function can make.
    'conv' denotes how the base10 values will be converted (should include padding).
    For a binary string, 'conv' should be 'paddedBinary'.
    """
    limb = ''

    if limbsPermitted <= 0:
        return limb

    chld = random.randint(0, min([MAX_CHILDREN-1, limbsPermitted-1]))
    angl = random.randint(0, MAX_ANGLE)
    leng = random.randint(0, MAX_LENGTH)
    freq = random.randint(0, MAX_FREQUENCY)
    ampl = random.randint(0, MAX_AMPLITUDE)
    phas = random.randint(0, MAX_PHASE)

    limb += conv(chld, 3) + \
            conv(angl, 5) + \
            conv(leng, 4) + \
            conv(freq, 4) + \
            conv(ampl, 4) + \
            conv(phas, 4)
    limbsPermitted -= 1

    for i in xrange(chld):
        childrenAfterThis = chld - i - 1
        chld_permitted = random.randint(1, limbsPermitted - childrenAfterThis)
        limbsPermitted -= chld_permitted
        limb += makeLimb(chld_permitted)
    
    return limb


def fitness_distanceMovedPerTime(humperdink):
    sim = subprocess.Popen(#'./simulator -i 1000 ',    # uncomment this line to see the organism
                           './simulator -i 1000 -g',
                           shell  = True, 
                           stdout = subprocess.PIPE, 
                           stdin  = subprocess.PIPE)

    coords = sim.communicate(humperdink.JSONObject)

    if coords == ('', None):
        return -100

    coords = coords[0].split('\n')[-2].split(',')
    x_coord = coords[0][1:]


    if x_coord == 'nan':
	return -100
    return float(x_coord)


def viewSimulation(humperdink):
    print 'sumiluaosnetuhaeo'
    sim = subprocess.Popen('./simulator -s 10',
                           shell = True,
                           stdin = subprocess.PIPE)
    sim.communicate(humperdink.JSONObject)



if __name__ == '__main__':
    main()

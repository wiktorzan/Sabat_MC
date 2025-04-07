from multiprocessing.dummy import Pool as PoolThread
from random import randint
from time import sleep
from os import system

def sim(waitTime):
    sleep(waitTime)
    system("./Sabat_exe macros/neutron14MeV.mac initConfig.dat >>../DATA/parallel_out" + str(waitTime))

processNumbersToSim = 20
numbers = set()
while len(numbers) < processNumbersToSim:
    numbers.add(randint(0, processNumbersToSim+1))

pool = PoolThread(6)
results = pool.map(sim, numbers)

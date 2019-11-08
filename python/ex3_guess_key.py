from random import random

def GetGuess():
    tmp = int(input("Please input your guess\n"))
    return tmp


key = int(random()*100)

print("guess the number in (0,100])")

guess = GetGuess()

while guess != key:
    if guess > key:
        print("Too large!")
    if guess < key:
        print("Too small!")
    guess = GetGuess()

print("Good job!")

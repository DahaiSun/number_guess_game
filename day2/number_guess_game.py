import random
random_number = random.randrange(1,21)

while True:
    guess = input("I have a number in mind between 1-20. What do you think it is? ")
    if guess < random_number:
        print("too small, try again") 
    if guess > random_number:
        print("too big, try again")
    if guess == random_number:
        print("you are right")
        
print(number)
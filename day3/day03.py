import random
random_number = random.randint(1,21)

guess_time= 0
while True:
    guess = int(input("I have a int number in mind between 1-20. What do you think it is? "))
    guess_time += 1
    if guess < random_number:
        print("too small, try again") 
    if guess > random_number:
        print("too big, try again")
    if guess == random_number:
        print("you are right")
        break
    
print(f"number of guesses: {guess_time}")

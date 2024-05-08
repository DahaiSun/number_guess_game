import random

def gernerate_random_number():
    return random.randint(1,21)

def obtain_user_input():
    return input("I have an int number in mind between 1-20. What do you think it is? (n-new game; s-cheat; x-exit)").lower()

def guess_process(guess, random_number):
    if guess < random_number:
        print("too small, try again")
    elif guess > random_number:
        print("too big, try again")
    elif guess == random_number:
        print("you are right")
        return True
    return False

def main():
    while True:
        random_number = gernerate_random_number()
        guess_time = 0
        while True:
            user_input = obtain_user_input()

            if user_input == "n":
                print("new game")
                break
                
            elif user_input == "s":
                print (f"the secret number is {random_number}")
                continue
            elif user_input == "x":
                return 
            
            try: 
                guess = int(user_input)
                guess_time += 1
                if guess_process(guess, random_number):
                    break
            except ValueError:
                print("Hi bro, that's not an integer")
        print( f"you tried {guess_time} times in the last game")
main()

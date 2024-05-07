import random

def main():

    while True:#big loop for the new game
        random_number = random.randint(1,21)
        guess_time = 0

        while True:#small loop for the guess          
            user_input = input("I have an int number in mind between 1-20. What do you think it is? (n-new game; s-cheat; x-exit)")
            
            if user_input.lower() == 'n':
                print("new game")
                break
            elif user_input.lower() == 's':
                print(f"the secret number is {random_number}")
                continue
            elif user_input.lower() == 'x':
                return#leave game()
            
            try:
                guess = int(user_input)
            except ValueError:
                print("Hi bro, that's not a integer")               
                    
            guess_time += 1
            if guess < random_number:
                print("too small, try again") 
            elif guess > random_number:
                print("too big, try again")
            elif guess == random_number:
                print("you are right")
                break
            
        print(f"number of guesses: {guess_time}")     

main()
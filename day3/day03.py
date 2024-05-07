import random

def main():

    def game(cheat=False):
        random_number = random.randint(1,21)
        if cheat:
            print(f"Hi loser, the number is: {random_number}")
        guess_time = 0
        while True:
            try:
                guess = int(input("I have an int number in mind between 1-20. What do you think it is? "))
            except ValueError:
                print("Hi bro, that is not a interger")
                continue
            
            
            guess_time += 1
            if guess < random_number:
                print("too small, try again") 
            elif guess > random_number:
                print("too big, try again")
            elif guess == random_number:
                print("you are right")
                break
            
        print(f"number of guesses: {guess_time}")
       

    game ()

    def new_game():
        while True:
            command = input("Play again? x-exit, n-new game, s-tell you the number")
            if command.lower() == "n":
                game()
            elif command.lower() == "s":
                game(cheat=True)
            elif command.lower() == "x":
                break
    new_game()
    
main()
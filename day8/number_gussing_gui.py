import tkinter as tk
from tkinter import messagebox
import random

def generate_random_number():
    return random.randint(1, 21)

def new_game():
    global random_number, guess_count
    random_number = generate_random_number()
    guess_count = 0
    guess_entry.delete(0, tk.END)
    feedback_label.config(text="Guess a number between 1 and 20")
    guess_count_label.config(text="Guess count: 0")

def show_number():
    top = tk.Toplevel(root)
    top.geometry("300x100")
    top.title("Secret Number")
    message = tk.Label(top, text=f"The secret number is {random_number}")
    message.pack(expand=True, fill=tk.BOTH, pady=20)

def check_guess():
    global guess_count
    try:
        guess = int(guess_entry.get())
        guess_count += 1
        guess_count_label.config(text=f"Guess count: {guess_count}")

        if guess < random_number:
            feedback_label.config(text="Too small, try again")
        elif guess > random_number:
            feedback_label.config(text="Too big, try again")
        else:
            feedback_label.config(text="You are right!")
    except ValueError:
        feedback_label.config(text="That's not an integer, try again")

def exit_game():
    root.destroy()

# Initialize the main window
root = tk.Tk()
root.title("Number Guessing Game")
root.geometry("300x300")
root.resizable(False, False)

# Create widgets
guess_entry = tk.Entry(root, width=20)
guess_entry.pack(pady=10)

submit_button = tk.Button(root, text="Submit Guess", command=check_guess)
submit_button.pack(pady=5)

new_game_button = tk.Button(root, text="New Game", command=new_game)
new_game_button.pack(pady=5)

show_number_button = tk.Button(root, text="Show Number", command=show_number)
show_number_button.pack(pady=5)

exit_button = tk.Button(root, text="Exit", command=exit_game)
exit_button.pack(pady=5)

feedback_label = tk.Label(root, text="Guess a number between 1 and 20")
feedback_label.pack(pady=10)

guess_count_label = tk.Label(root, text="Guess count: 0")
guess_count_label.pack(pady=5)

# Start a new game initially
new_game()

# Start the main loop
root.mainloop()

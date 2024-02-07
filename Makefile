all:
	python3 ./generator_sequences/generator.py
	g++ -std=c++20 -Wall -Wextra -Werror ./src/*.cpp -o ./olp_asm
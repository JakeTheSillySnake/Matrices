.PHONY: s21_matrix.a test gcov_report check rebuild clean
CC = gcc -Wall -Wextra -Werror -std=c11 $(shell pkg-config --cflags check)
GCOVFLAGS = -fprofile-arcs -ftest-coverage
OS=$(shell uname -s)
LIBS= $(shell pkg-config --libs check)
LINUX = -pthread -lcheck_pic -pthread -lrt -lm -lsubunit 
HEADERS = s21_matrix.h matrix_extra.h
SOURCES = $(wildcard *.c)
SRC= s21_matrix.c matrix_extra.c

OBJ_LIBRARY := $(patsubst %.c, %.o, $(SRC))
ifeq ($(OS), Linux)
 CC += -D OS_LINUX
endif
ifeq ($(OS), Darwin)
 CC += -D OS_MAC
endif

%.o: %.c $(HEADER)
	$(CC) -c $< -o $@

all: clean s21_matrix.a test gcov_report

s21_matrix.a: $(OBJ_LIBRARY) $(HEADERS)
	ar rcs s21_matrix.a $(OBJ_LIBRARY)
	ranlib s21_matrix.a
	rm -rf *.o

test: s21_matrix.a tests.c
ifeq ($(OS), Darwin)
	$(CC) tests.c s21_matrix.a -o test.out $(LIBS) -lm
else
	$(CC) tests.c s21_matrix.a -o test.out $(LIBS) $(LINUX)
endif
	./test.out

gcov_report: clean_test
ifeq ($(OS), Darwin)
	$(CC) $(GCOVFLAGS) tests.c ${SRC} -o test.out $(LIBS) -lm
else
	$(CC) $(GCOVFLAGS) tests.c ${SRC} -o test.out $(LIBS) $(LINUX)
endif
	./test.out
	lcov -t "stest" -o s21_test.info -c -d .
	genhtml -o report s21_test.info
ifeq ($(OS), Linux)
	links ./report/index.html
else
	open -a safari ./report/index.html
endif

check:
	cppcheck --enable=all --suppress=missingIncludeSystem ./
	clang-format -i *.c
	clang-format -i *.h

leaks:
ifeq ($(OS), Darwin)
	$(CC) $(SRC) $(HEADERS) tests.c -o main.out $(LIBS) -lm
	leaks --atExit -- ./main.out
else
	$(CC) $(SRC) $(HEADERS) tests.c -o main.out $(LIBS) $(LINUX)
	valgrind --tool=memcheck --leak-check=yes ./main.out
endif
	rm -f main.out

rebuild: clean all

clean:
	rm -f *.o *.a *.gcno *.gcda *.info *.out *.gcov
	rm -rf test
	rm -rf ./report
	rm -rf *.dSYM

clean_test:
	rm -rf *.gcda
	rm -rf *.gcno
	rm -rf *.info
	rm -rf test
	rm -rf report

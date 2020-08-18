OBJ = \
source/needleman_wunsch.o \
source/needleman_wunsch_affine_gaps.o \
source/smith_waterman.o \
source/string_functions.o \
main/align_pairs.o

COMPILER = g++

NAME = alignme1.0.exe

CPPFLAGS = -Wall -O3 -fmessage-length=0 -Wno-deprecated 

.o:
	$(COMPILER) -c $(CPPFLAGS) -o $@  $< 

all: $(NAME)

$(NAME): $(OBJ)
	$(COMPILER) -o $(NAME) $(OBJ)

clean:
	rm -f $(NAME)
	rm -f $(OBJ)


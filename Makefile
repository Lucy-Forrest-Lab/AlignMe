OBJ = \
./source/needleman_wunsch.o \
./source/needleman_wunsch_affine_gaps.o \
./source/reader.o \
./source/score_profile_similarity.o \
./source/score_profile_similarity_linear_normalized.o \
./source/score_sequence_similarity.o \
./source/score_sequence_similarity_profile_dependent.o \
./source/sequence.o \
./source/smith_waterman.o \
./source/string_functions.o \
./main/align_pairs.o

CXX = g++

NAME = alignme1.2.exe

CXXFLAGS = -Wall -O3 -fmessage-length=0 -Wno-deprecated

LDLIBS = -L/usr/local/lib
LDFLAGS = -static

$(NAME): $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -o $(NAME) $(OBJ)

clean:
	rm -f $(NAME)
	rm -f $(OBJ)


compiler		=mpicxx

ifeq ($(mode),debug)
	cflags		+=-std=c++17 -ggdb -Wno-unused-result -fsanitize=address -fno-omit-frame-pointer -march=native -fopenmp -I ./include
	linkflags   += -fsanitize=address -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -fopenmp
else
	cflags		+=-std=c++17 -Ofast	-I ./include
	linkflags	+=
endif

SrcDir			=./src
HeadDir			+=./include
ObjDir			=./obj
OutDir			=./output/results ./output/time_record/compute ./output/time_record/total ./output/load ./output/grid
source			=$(foreach dir,$(SrcDir),$(wildcard $(dir)/*.cpp))
head			=$(foreach dir,$(HeadDir),$(wildcard $(dir)/*.h))
object			=$(patsubst %.cpp,$(ObjDir)/%.o,$(notdir $(source)))
output			=$(foreach dir,$(OutDir),$(wildcard $(dir)/ *.txt $(dir)/*.dat))
target 			=EU3D
NO_COLOR		=\033[0m
OK_COLOR		=\033[32,01m

$(target):$(object) $(head)
	$(compiler) -o $(target) $(object) $(linkflags)
	@printf "$(OK_COLOR)Compiling Is Successful!\nExecutable File: $(target) $(NO_COLOR)\n"

$(ObjDir)/%.o:$(SrcDir)/%.cpp $(head)
	$(compiler) -c $(cflags) $< -o $@

.PHONY:clean,cleano
clean:
	-rm $(object) $(target)
cleano:
	-rm $(output)
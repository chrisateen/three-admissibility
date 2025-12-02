# 3-admissibility
Rust implementation of the 3-admissibility algorithm

## Link to paper
https://arxiv.org/abs/2512.01121

## To run the program
1. This program can be run using the collection of networks found in the network-corpus repo
2. Download the network-corpus repo and run the following command to run the program on a specific network
`three-admissibility <NAME_OF_NETWORK> <P_VALUE_TO_START_THE_SEARCH> <DIR_TO_THE_NETWORK>`

For example
`three-admissibility windsurfers 11 ../network-corpus/networks`

### Saving ordering to a file
The 3-admissibility ordering for a graph can be saved to a txt.gz file by using the save command.

`three-admissibility <NAME_OF_NETWORK> <P_VALUE_TO_START_THE_SEARCH> <DIR_TO_THE_NETWORK> save <OPTIONAL_DIR_WHERE_TO_SAVE_THE_FILE>`

For example the below would save orderings for the graph in /Users/username/orderings/windsurfers.txt.gz
`three-admissibility windsurfers 11 ../network-corpus/networks /Users/username/orderings`

If no directory is included a results directory will be created in the directory where the program is running and save the results there

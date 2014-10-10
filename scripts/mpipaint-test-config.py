TilePadding = 1

M = [
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1]]

PIXEL_WIDTH = 2200
CHUNK_SIZE = 1024 * 1024 * 12
SMALL_IMAGE = False
POS_BLOCK = "0/Position"
SML_BLOCK = "0/SmoothingLength"
DATA_BLOCK =  ["0/Mass", 
        "0/StarFormationRate",
        (ieye2Tmass, "0/InternalEnergy", "0/ElectronAbundance", "0/Mass")
]

#NEAR=0
#FAR=1000

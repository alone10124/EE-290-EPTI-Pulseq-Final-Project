# EE 290 Final Project - EPTI Pulse Sequence Improvements

## Folder Structure: 

./
├── config/                # .yaml (or .json) configuration files that will determine flow 
├── out/                   # Plots, results of simulations            
├── sequences/             # .seq files. Can be imported to GE scanner or KomaMRI sim scripts
│   └── scripts/           # PyPulseq code that will generate .seq files
├── sim/                   # Simulation scripts, using KomaMRI Python wrapper
├── .gitignore             # Files ignored by git
├── LICENSE                # Project license details
└── README.md              # Project overview (this file)
# EE 290 Final Project - EPTI Pulse Sequence Improvements

## Folder Structure: 

./
├── config/                # .yaml (or .json) configuration files that will determine flow 
├── out/                   # Plots, results of simulations
├── seq_gen/               # PyPulseq code that will generate .seq files              
├── seqs/                  # .seq files generated from seq_gen/. Can be imported to GE scanner or KomaMRI sim scripts                  
├── sim/                   # Simulation scripts, using KomaMRI Python wrapper
├── .gitignore             # Files ignored by git
├── LICENSE                # Project license details
└── README.md              # Project overview (this file)
# Data Analysis Project

## Overview

This repository contains data files, configuration files, and documentation related to an indentation analysis project. The project appears to involve experimental data collection and analysis of indentation parameters including radius, size, and velocity measurements.

## Project Structure

### Configuration Files

- **`Configuration.xls`** - Main configuration file containing experimental setup parameters
- **`Configuration-Iron.xls`** - Configuration file specific to iron material testing

### Data Files

- **`Radius.xlsx`** - Contains radius measurement data
- **`Size.xlsx`** - Contains size/dimension measurement data
- **`Velocity.xlsx`** - Contains velocity measurement data

### Results

- **`GoodResult.xlsx`** - Processed results and analysis output

### Scripts

- **`plot_indentation.py`** - Python script to plot indentation data with dual y-axes:
  - Left axis: Indentation time (Timestep*10000)[ps] vs Force Z(*1.602) [nN]
  - Right axis: Indentation time (Timestep*10000)[ps] vs Hardness [GPa]


*Note: Temporary files (starting with ~$) are automatically created by Microsoft Office applications when files are open and can be safely ignored.*

## Usage

### Installation

To use the plotting script, first install the required Python packages:

```bash
pip install -r requirements.txt
```

### Plotting Indentation Data

To generate plots from `GoodResult.xlsx`:

```bash
python plot_indentation.py
```

This will create a dual-axis plot showing:
- **Left axis**: Force Z (*1.602) [nN] vs Indentation Time (Timestep*10000) [ps]
- **Right axis**: Hardness [GPa] vs Indentation Time (Timestep*10000) [ps]

The plot will be saved as `indentation_plot.png` and displayed on screen.

### Other Usage

1. **Configuration Setup**: Review and modify configuration files (`Configuration.xls` or `Configuration-Iron.xls`) to set up experimental parameters.

2. **Data Analysis**: 
   - Load data from `Radius.xlsx`, `Size.xlsx`, and `Velocity.xlsx`
   - Process and analyze the data
   - Generate results in `GoodResult.xlsx`

3. **Documentation Review**: 
   - Refer to `indentation_NG_TL.pdf` for technical details
   - Review `Paper1_Review_indentation.pptx` for presentation materials

## File Formats

- **Excel Files (.xls, .xlsx)**: Data and configuration files
- **Python Scripts (.py)**: Data analysis and plotting scripts
- **PDF (.pdf)**: Technical documentation
- **PowerPoint (.pptx)**: Presentation materials
- **Image Files (.png)**: Generated plots and visualizations

## Notes

- Ensure Microsoft Excel or compatible software is installed to open `.xls` and `.xlsx` files
- Temporary lock files (starting with `~$`) are created automatically and should not be manually edited
- Always close Excel/PowerPoint files before committing changes to avoid leaving temporary files

## Contact

For questions or issues related to this project, please refer to the documentation files or contact the project maintainer.


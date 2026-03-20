# Contents: /user_definitions

This folder contains user-specific configuration data in `.csv` format. The file is automatically created during the first initialization of the repository in MATLAB.

## Fields

- **user** → Your Pulseq username  
- **lab** → Your MRI lab or institutional affiliation  
- **scanner** → Your default MRI scanner  
- **path** → Backup path where important metadata is stored (e.g., in the `PULSEQ` struct)

This configuration ensures reproducibility and consistency across different machines or users in a lab environment.

_Maximilian Gram: 14.07.2025_
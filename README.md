# Flash_Dock

To use Flash_Dock, you need to download the `unimol_docking_v2_240517.pt` model file.

## Requirements

- Use uv

```bash
uv sync
```

- Use pip

```bash
pip install -r requirements.txt
```

## Steps:

1. Download the model weights and install the required dependencies from the following link:  
   [https://github.com/deepmodeling/Uni-Mol/tree/main/unimol_docking_v2](https://github.com/deepmodeling/Uni-Mol/tree/main/unimol_docking_v2)

2. Place the file in the directory:  
   `./others/Uni-Mol/unimol_docking_v2/`

3. Start the application using the following command:  

   ```bash
   streamlit run FlashDock.py
   ```

   - Use uv:

      ```bash
      uv run streamlit run FlashDock.py
      ```

Ensure the file is in the correct location for the program to function properly.

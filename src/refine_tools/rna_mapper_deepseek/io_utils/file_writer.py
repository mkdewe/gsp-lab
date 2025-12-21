"""Zapis wyników do plików"""

from utils.logger import Logger

class FileWriter:
    @staticmethod
    def save_mapping(filename, mapping, metadata=None):
        with open(filename, 'w') as f:
            if metadata:
                f.write(f"# Target: {metadata.get('target', '')}\n")
                f.write(f"# Model: {metadata.get('model', '')}\n\n")
            for model_id, target_id in sorted(mapping.items()):
                f.write(f"{model_id[0]},{model_id[1]} -> {target_id[0]},{target_id[1]}\n")
        Logger.info(f"Zapisano: {filename}")
    
    @staticmethod
    def save_missing_atoms(filename, missing_atoms):
        if not missing_atoms: 
            return
        with open(filename, 'w') as f:
            for report in missing_atoms:
                f.write(f"{report['model_chain']},{report['model_residue']} -> "
                       f"{report['target_chain']},{report['target_residue']}: "
                       f"{','.join(report['missing_atoms'])}\n")
        Logger.info(f"Zapisano: {filename}")
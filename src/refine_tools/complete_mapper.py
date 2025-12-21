import numpy as np
import warnings
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.Align import PairwiseAligner
from collections import defaultdict
from itertools import permutations
import argparse
import os
import sys

# WYŁĄCZENIE WSZYSTKICH WARNINGÓW
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=PDBConstructionWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore")

class CompleteRNAMapper:
    def __init__(self):
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -0.5
        self.aligner.extend_gap_score = -0.1
        
        # Sprawdź czy scipy jest dostępny
        try:
            from scipy.optimize import linear_sum_assignment
            self.has_scipy = True
        except ImportError:
            self.has_scipy = False
            print("UWAGA: scipy nie jest zainstalowany. Mapowanie geometryczne nie będzie dostępne.")
    
    def main():
        parser = argparse.ArgumentParser(description='Mapowanie RNA między targetem a modelem')
        parser.add_argument('target_pdb', help='Ścieżka do pliku PDB target')
        parser.add_argument('model_pdb', help='Ścieżka do pliku PDB modelu')
        parser.add_argument('--output', required=True, help='Ścieżka do wyjściowego pliku PDB')
        parser.add_argument('--remove-extra-atoms', action='store_true', 
                        help='Usuń nadmiarowe atomy z modelu')
        
        args = parser.parse_args()
        
        mapper = CompleteRNAMapper()
        result = mapper.process_all(
            args.target_pdb,
            args.model_pdb,
            args.output,
            remove_extra_atoms=args.remove_extra_atoms
        )
        
        # Sprawdź czy utworzono solution
        residue_mapping, missing_atoms_report, rmsd, changed_count, removed_atoms_count, removed_residues_count, atom_comparison_log, total_extra_atoms, total_missing_atoms, solution_output = result
        if solution_output and os.path.exists(solution_output):
            print(f"Utworzono plik solution: {solution_output}")
        else:
            print("Nie utworzono pliku solution")

    if __name__ == "__main__":
        main()


#region Proste

    def count_rna_residues(self, structure):
        """Liczy liczbę reszt RNA w strukturze."""
        count = 0
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == ' ':  # tylko standardowe reszty
                        resname = residue.resname.strip()
                        if resname in ['A','C','G','U']:
                            count += 1
        return count

    def extract_chain_sequences(self, structure):
        """Ekstrahuje sekwencje dla każdego łańcucha."""
        chain_sequences = {}
        chain_residues = {}
        
        for model in structure:
            for chain in model:
                seq = ""
                residues = []
                for residue in chain:
                    if residue.id[0] == ' ':  # tylko standardowe reszty
                        resname = residue.resname.strip()
                        if resname in ['A','C','G','U']:
                            seq += resname
                            residues.append(residue)
                if seq:  # tylko niepuste sekwencje
                    chain_sequences[chain.id] = seq
                    chain_residues[chain.id] = residues
        
        return chain_sequences, chain_residues
    
    def compute_chain_centroids(self, structure, chain_residues):
        """Oblicza centroidy (środki geometryczne) dla każdego łańcucha."""
        chain_centroids = {}
        
        for chain_id, residues in chain_residues.items():
            if not residues:
                continue
                
            coords = []
            for residue in residues:
                # Używamy atomu P lub C4' jako reprezentanta reszty
                atom = None
                try:
                    atom = residue['P']
                except KeyError:
                    try:
                        atom = residue['C4\'']
                    except KeyError:
                        # Jeśli brak obu atomów, spróbuj dowolnego atomu
                        for atom_obj in residue:
                            if atom_obj.name not in ['H', 'OP1', 'OP2']:  # Unikaj atomów wodoru i tlenów fosforanowych
                                atom = atom_obj
                                break
                
                if atom is not None:
                    coords.append(atom.get_coord())
            
            if coords:
                centroid = np.mean(coords, axis=0)
                chain_centroids[chain_id] = centroid
        
        return chain_centroids
    
    def get_target_chain_order(self, target_structure):
        """Zwraca kolejność łańcuchów z targetu."""
        chain_order = []
        for model in target_structure:
            for chain in model:
                chain_order.append(chain.id)
        return chain_order
#endregion 
    
    def find_optimal_chain_mapping_by_geometry(self, target_centroids, model_centroids, target_sequences, model_sequences):
        """Znajduje optymalne mapowanie łańcuchów używając kombinacji sekwencji i geometrii."""
        from scipy.optimize import linear_sum_assignment
        
        # Jeśli liczba łańcuchów się nie zgadza, używamy tylko mapowania sekwencyjnego
        if len(target_centroids) != len(model_centroids):
            print("  Liczba łańcuchów się nie zgadza, używam mapowania sekwencyjnego")
            return self.find_chain_mapping_by_sequence(target_sequences, model_sequences)
        
        target_chains = sorted(list(target_centroids.keys()))
        model_chains = sorted(list(model_centroids.keys()))
        
        # 1. Macierz tylko geometryczna (ignoruje sekwencje)
        geo_matrix = np.zeros((len(target_chains), len(model_chains)))
        for i, t_chain in enumerate(target_chains):
            for j, m_chain in enumerate(model_chains):
                geo_matrix[i, j] = np.linalg.norm(
                    target_centroids[t_chain] - model_centroids[m_chain]
                )
        
        # 2. Macierz kombinowana (sekwencja + geometria)
        seq_matrix = np.zeros((len(target_chains), len(model_chains)))
        for i, t_chain in enumerate(target_chains):
            for j, m_chain in enumerate(model_chains):
                if t_chain in target_sequences and m_chain in model_sequences:
                    seq_matrix[i, j] = self.calculate_sequence_similarity(
                        target_sequences[t_chain], model_sequences[m_chain]
                    )
        
        # Normalizacja
        if np.max(geo_matrix) > 0:
            geo_matrix_normalized = geo_matrix / np.max(geo_matrix)
        else:
            geo_matrix_normalized = geo_matrix
        
        # Używamy ODWROTNOŚCI podobieństwa sekwencji jako kary
        penalty_matrix = 1.0 - seq_matrix
        
        # Połącz: geometria + kara za różne sekwencje
        sequence_penalty_weight = 0.3
        combined_matrix = geo_matrix_normalized + (sequence_penalty_weight * penalty_matrix)
        
        # Znajdź minimalny koszt (najbliższe centroidy + najmniejsza kara sekwencyjna)
        row_ind_geo, col_ind_geo = linear_sum_assignment(geo_matrix)
        row_ind_combined, col_ind_combined = linear_sum_assignment(combined_matrix)
        
        # Oblicz całkowity koszt dla obu rozwiązań
        geo_cost = geo_matrix[row_ind_geo, col_ind_geo].sum()
        combined_cost = combined_matrix[row_ind_combined, col_ind_combined].sum()
        
        print(f"  Porównanie kosztów:")
        print(f"    Tylko geometria: {geo_cost:.3f}")
        print(f"    Geometria + sekwencja: {combined_cost:.3f}")
        
        # Wybierz rozwiązanie z niższym kosztem
        if geo_cost <= combined_cost * 1.1:  # Tolerancja 10%
            print("  Wybrano mapowanie oparte TYLKO na geometrii")
            row_ind, col_ind = row_ind_geo, col_ind_geo
            matrix_used = "geometry"
        else:
            print("  Wybrano mapowanie kombinowane (geometria + sekwencja)")
            row_ind, col_ind = row_ind_combined, col_ind_combined
            matrix_used = "combined"
        
        # Tworzymy mapowanie
        chain_mapping = {}
        for i, j in zip(row_ind, col_ind):
            target_chain = target_chains[i]
            model_chain = model_chains[j]
            chain_mapping[target_chain] = model_chain
            
            # Oblicz metryki dla informacji
            geo_dist = geo_matrix[i, j]
            seq_sim = seq_matrix[i, j] if i < seq_matrix.shape[0] and j < seq_matrix.shape[1] else 0
            
            print(f"    {target_chain} -> {model_chain} (odległość: {geo_dist:.2f}Å, podob. sekwencji: {seq_sim:.3f})")
        
        print(f"  Użyta macierz: {matrix_used}")
        
        return chain_mapping
    
    def find_chain_mapping_by_sequence(self, target_sequences, model_sequences):
        """Oryginalne mapowanie oparte tylko na sekwencji - BEZ arbitralnych progów."""
        chain_mapping = {}
        used_model_chains = set()
        
        print("  Mapowanie oparte na sekwencji:")
        
        for target_chain, target_seq in target_sequences.items():
            best_score = -float('inf')
            best_model_chain = None
            
            for model_chain, model_seq in model_sequences.items():
                if model_chain in used_model_chains:
                    continue
                
                score = self.calculate_sequence_similarity(target_seq, model_seq)
                
                # ZAWSZE wybierz najlepszy score, nawet jeśli niski!
                if score > best_score:
                    best_score = score
                    best_model_chain = model_chain
            
            if best_model_chain:
                chain_mapping[target_chain] = best_model_chain
                used_model_chains.add(best_model_chain)
                print(f"    {target_chain} -> {best_model_chain} (score: {best_score:.3f})")
            else:
                print(f"    OSTRZEŻENIE: Nie znaleziono mapowania dla łańcucha {target_chain}")
        
        return chain_mapping
    
    def find_chain_mapping(self, target_sequences, model_sequences):
        """
        Znajduje optymalne mapowanie łańcuchów między targetem a modelem,
        najpierw próbując dopasować przestrzennie, potem sekwencyjnie.
        """
        chain_mapping = {}
        used_model_chains = set()
        print(" Szczegółowe porównanie sekwencji + geometrii:")

        # POPRAWIONE: sprawdź czy atrybuty istnieją i czy są niepuste
        use_spatial = (
            hasattr(self, '_last_target_residues') and 
            hasattr(self, '_last_model_residues') and
            self._last_target_residues is not None and
            self._last_model_residues is not None and
            len(self._last_target_residues) > 0 and
            len(self._last_model_residues) > 0
        )
        
        print(f"  DEBUG: use_spatial = {use_spatial}")
        if use_spatial:
            print(f"  DEBUG: target_residues keys: {list(self._last_target_residues.keys())}")
            print(f"  DEBUG: model_residues keys: {list(self._last_model_residues.keys())}")
        
        # Najpierw szybki pass: dla każdego target_chain wybierz najlepszy model_chain wg kombinacji spatial i sekwencji.
        for target_chain, target_seq in target_sequences.items():
            best_score = -1e9
            best_model_chain = None
            
            for model_chain, model_seq in model_sequences.items():
                if model_chain in used_model_chains:
                    continue

                # sekwencyjny score (0..1)
                seq_score = self.calculate_sequence_similarity(target_seq, model_seq)

                # spatial score (mniejszy RMSD = lepszy) -> zamienimy na skalowany score
                spatial_score = 0.0
                if use_spatial:
                    t_res = self._last_target_residues.get(target_chain, [])
                    m_res = self._last_model_residues.get(model_chain, [])
                    print(f"  DEBUG: target_chain={target_chain}, t_res length={len(t_res)}")
                    print(f"  DEBUG: model_chain={model_chain}, m_res length={len(m_res)}")
                    rmsd = self.chain_spatial_similarity(t_res, m_res)
                    if rmsd != float('inf'):
                        # przekształcamy rmsd na score: mniejszy rmsd -> wyższy score
                        spatial_score = 1.0 / (1.0 + rmsd)  # zakres (0,1]
                    else:
                        spatial_score = 0.0

                # Połączone kryterium: dajemy większy ciężar geometrii (np. 0.7) i sekwencji (0.3)
                combined_score = 0.0
                if use_spatial:
                    combined_score = 0.7 * spatial_score + 0.3 * seq_score
                else:
                    combined_score = seq_score

                print(f" {target_chain}(target) vs {model_chain}(model): seq={seq_score:.3f} spatial={'{:.3f}'.format(spatial_score) if use_spatial else 'N/A'} combined={combined_score:.3f}")

                # Wybierz model_chain jeśli score lepszy i spełnia minimalny próg
                # ZMIANA: obniż próg do 0.1 dla testów
                if combined_score > best_score and combined_score > 0.1:
                    best_score = combined_score
                    best_model_chain = model_chain

            if best_model_chain:
                chain_mapping[target_chain] = best_model_chain
                used_model_chains.add(best_model_chain)
                print(f" Mapowanie: {target_chain} -> {best_model_chain} (score: {best_score:.3f})")
            else:
                print(f"  OSTRZEŻENIE: Nie znaleziono mapowania dla łańcucha {target_chain}")

        # Jeśli nie znaleziono mapowania (np. gdy sekwencje są bardzo różne) spróbuj permutacji wg RMSD (jeśli dostępne reszty)
        if (not chain_mapping or len(chain_mapping) < len(target_sequences)) and use_spatial:
            print("  Próba mapowania permutacyjnego...")
            perm_map, perm_rmsd = self.find_permutation_by_rmsd(self._last_target_residues, self._last_model_residues)
            if perm_map:
                # perm_map: target_chain -> model_chain
                # zachowaj tylko nowe pary które nie kolidują z istniejącym mappingiem
                for t, m in perm_map.items():
                    if t not in chain_mapping and m not in used_model_chains:
                        chain_mapping[t] = m
                        used_model_chains.add(m)
                print(f"Zastosowano permutacyjne mapowanie wg RMSD (sum RMSD={perm_rmsd:.3f}): {perm_map}")

        # Jeśli nadal brak mapowania, użyj tylko sekwencji
        if not chain_mapping:
            print("  Mapowanie geometryczne nie dało wyników, używam mapowania tylko po sekwencji")
            chain_mapping = self.find_chain_mapping_by_sequence(target_sequences, model_sequences)
        
        return chain_mapping

    def calculate_sequence_similarity(self, seq1, seq2):
        """Oblicza podobieństwo między sekwencjami."""
        if seq1 == seq2:
            return 1.0
        
        try:
            alignments = self.aligner.align(seq1, seq2)
            if alignments:
                best_alignment = alignments[0]
                alignment_length = len(best_alignment[0])
                matches = sum(1 for a, b in zip(best_alignment[0], best_alignment[1]) if a == b)
                return matches / alignment_length
        except Exception:
            pass
        
        return 0.0

    def chain_spatial_similarity(self, residues1, residues2):
        """
        Liczy RMSD przestrzenny między dwoma łańcuchami używając atomów 'P' (fallback: C4').
        Zwraca float (RMSD) lub float('inf') jeśli za mało punktów.
        """
        from Bio.PDB import Superimposer

        atoms1 = []
        atoms2 = []

        # Dopasowujemy pary według kolejności reszt (zip krótszy)
        for r1, r2 in zip(residues1, residues2):
            # Preferuj atom 'P', jeśli brak to 'C4''
            a1 = None
            a2 = None
            try:
                a1 = r1['P']
            except KeyError:
                try:
                    a1 = r1["C4'"]
                except KeyError:
                    # Spróbuj znaleźć jakikolwiek atom
                    for atom in r1:
                        if atom.name not in ['H', 'OP1', 'OP2']:
                            a1 = atom
                            break
            try:
                a2 = r2['P']
            except KeyError:
                try:
                    a2 = r2["C4'"]
                except KeyError:
                    # Spróbuj znaleźć jakikolwiek atom
                    for atom in r2:
                        if atom.name not in ['H', 'OP1', 'OP2']:
                            a2 = atom
                            break

            if a1 is not None and a2 is not None:
                atoms1.append(a1)
                atoms2.append(a2)

        print(f"    DEBUG chain_spatial_similarity: znaleziono {len(atoms1)} par atomów")
        
        if len(atoms1) < 3:
            print(f"    DEBUG: Za mało atomów do obliczenia RMSD ({len(atoms1)} < 3)")
            return float('inf')

        sup = Superimposer()
        try:
            sup.set_atoms(atoms1, atoms2)
            print(f"    DEBUG: RMSD = {sup.rms:.3f}")
            return sup.rms
        except Exception as e:
            print(f"    DEBUG: Błąd superimposera: {e}")
            return float('inf')

    def find_permutation_by_rmsd(self, target_residues, model_residues, max_perm=8):
        """
        Szuka najlepszej permutacji modelowych łańcuchów dopasowujących się do łańcuchów targetu
        przez minimalizację sumy RMSD. Zwraca (chain_mapping_dict, best_total_rmsd).
        Jeśli liczba łańcuchów modelu > max_perm, zwraca ({}, float('inf')).
        """
        target_chains = list(target_residues.keys())
        model_chains = list(model_residues.keys())

        if len(model_chains) != len(target_chains):
            return {}, float('inf')

        if len(model_chains) > max_perm:
            # Nie robimy pełnych permutacji dla zbyt wielu łańcuchów
            return {}, float('inf')

        best_perm = None
        best_total = float('inf')

        for perm in permutations(model_chains):
            total = 0.0
            valid = True
            for t_chain, m_chain in zip(target_chains, perm):
                rmsd = self.chain_spatial_similarity(target_residues[t_chain], model_residues[m_chain])
                if rmsd == float('inf'):
                    valid = False
                    break
                total += rmsd
            if valid and total < best_total:
                best_total = total
                best_perm = perm

        if best_perm is None:
            return {}, float('inf')

        mapping = {}
        for t_chain, m_chain in zip(target_chains, best_perm):
            mapping[t_chain] = m_chain
        return mapping, best_total

    def optimize_chain_order_for_single_chain_model(self, target_sequences, model_sequence, target_residues, model_residues):
        """Optymalizuje kolejność nici w modelu jednonicowym dla lepszego dopasowania przestrzennego."""
        print("  Optymalizacja kolejności nici w modelu jednonicowym...")
        
        from itertools import permutations
        
        target_chains = list(target_sequences.keys())
        
        # Jeśli za dużo łańcuchów, użyj heurystyki zamiast wszystkich permutacji
        if len(target_chains) > 4:
            print(f"  Zbyt wiele łańcuchów ({len(target_chains)}), używanie heurystyki...")
            return self._heuristic_chain_ordering(target_sequences, model_sequence)
        
        best_order = None
        best_score = float('-inf')
        
        # Testuj różne kolejności łańcuchów targetu
        for order in permutations(target_chains):
            # Połącz sekwencje targetu w testowanej kolejności
            combined_target_seq = ''.join(target_sequences[chain] for chain in order)
            
            # Wyrównaj z sekwencją modelu
            try:
                alignments = self.aligner.align(combined_target_seq, model_sequence)
                if alignments:
                    aligned_target, aligned_model = alignments[0][0], alignments[0][1]
                    alignment_score = alignments[0].score
                    
                    # Dodaj kary za luki
                    gap_penalty = aligned_target.count('-') * 2 + aligned_model.count('-') * 2
                    final_score = alignment_score - gap_penalty
                    
                    if final_score > best_score:
                        best_score = final_score
                        best_order = order
            except Exception as e:
                continue
        
        print(f"  Najlepsza kolejność: {best_order} (score: {best_score:.2f})")
        return list(best_order) if best_order else target_chains

    def _heuristic_chain_ordering(self, target_sequences, model_sequence):
        """Heurystyczne ustalanie kolejności łańcuchów dla przypadków z wieloma łańcuchami."""
        print("  Używanie heurystycznej kolejności łańcuchów...")
        
        # Posortuj łańcuchy według długości (najdłuższy pierwszy)
        sorted_chains = sorted(target_sequences.keys(), 
                              key=lambda x: len(target_sequences[x]), 
                              reverse=True)
        
        # Spróbuj dopasować najdłuższy łańcuch do początku sekwencji modelu
        best_order = sorted_chains
        best_score = self.calculate_sequence_similarity(
            ''.join(target_sequences[chain] for chain in sorted_chains),
            model_sequence
        )
        
        # Testuj kilka wariantów
        test_orders = [sorted_chains]
        
        # Odwrócona kolejność
        test_orders.append(sorted_chains[::-1])
        
        # Próbuj umieścić krótkie łańcuchy na początku
        if len(sorted_chains) > 2:
            short_first = sorted(target_sequences.keys(), 
                               key=lambda x: len(target_sequences[x]))
            test_orders.append(short_first)
            test_orders.append(short_first[::-1])
        
        for order in test_orders:
            combined_seq = ''.join(target_sequences[chain] for chain in order)
            score = self.calculate_sequence_similarity(combined_seq, model_sequence)
            if score > best_score:
                best_score = score
                best_order = order
        
        print(f"  Heurystyczna kolejność: {best_order} (score: {best_score:.3f})")
        return best_order

    def create_residue_mapping_for_combined_chain(self, target_sequences, target_residues, model_chain, model_residues):
        """Tworzy mapowanie reszt gdy model ma jeden łańcuch, a target wiele - Z OPTYMALIZACJĄ KOLEJNOŚCI."""
        residue_mapping = {}
        missing_atoms_report = []
        
        print(f"  Mapowanie połączonego łańcucha {model_chain} na wiele łańcuchów targetu")
        
        # OPTYMALIZACJA: Znajdź najlepszą kolejność łańcuchów targetu
        model_seq = ''.join([r.resname.strip() for r in model_residues[model_chain]])
        optimal_chain_order = self.optimize_chain_order_for_single_chain_model(
            target_sequences, model_seq, target_residues, model_residues
        )
        
        # Połącz sekwencje i reszty targetu w OPTYMALNEJ kolejności
        combined_target_seq = ""
        target_residue_list = []
        
        for chain_id in optimal_chain_order:
            if chain_id in target_sequences:
                seq = target_sequences[chain_id]
                residues = target_residues[chain_id]
                for res in residues:
                    combined_target_seq += res.resname.strip()
                    target_residue_list.append((chain_id, res))

        print(f"    Optymalna kolejność łańcuchów: {optimal_chain_order}")
        print(f"    Połączony target: {combined_target_seq}")
        print(f"    Model: {model_seq}")

        # Reszta metody pozostaje bez zmian (alignowanie i mapowanie)
        try:
            alignments = self.aligner.align(combined_target_seq, model_seq)
            if alignments:
                best_alignment = alignments[0]
                aligned_target, aligned_model = best_alignment[0], best_alignment[1]

                print(f"    Wyrównanie target: {aligned_target}")
                print(f"    Wyrównanie model:  {aligned_model}")

                # Mapowanie po alignmencie
                target_idx = 0
                model_idx = 0

                for t_char, m_char in zip(aligned_target, aligned_model):
                    if t_char != '-' and m_char != '-':
                        if target_idx < len(target_residue_list) and model_idx < len(model_residues[model_chain]):
                            target_chain_id, target_res = target_residue_list[target_idx]
                            model_res = model_residues[model_chain][model_idx]

                            missing_atoms = self.check_missing_atoms(target_res, model_res)
                            if missing_atoms:
                                missing_atoms_report.append({
                                    'model_chain': model_chain,
                                    'model_residue': model_res.id[1],
                                    'target_chain': target_chain_id,
                                    'target_residue': target_res.id[1],
                                    'missing_atoms': missing_atoms
                                })

                            residue_mapping[(model_chain, model_res.id[1])] = (target_chain_id, target_res.id[1])
                            print(f"   {model_chain}{model_res.id[1]} -> {target_chain_id}{target_res.id[1]}")

                            target_idx += 1
                            model_idx += 1
                    elif t_char == '-':
                        if model_idx < len(model_residues[model_chain]):
                            print(f"    - Gap w targetcie, pomijam {model_chain}{model_residues[model_chain][model_idx].id[1]}")
                            model_idx += 1
                    elif m_char == '-':
                        if target_idx < len(target_residue_list):
                            print(f"    - Gap w modelu, pomijam {target_residue_list[target_idx][0]}{target_residue_list[target_idx][1].id[1]}")
                            target_idx += 1

                print(f"    Zmapowano {len(residue_mapping)} reszt w trybie połączonym")
            else:
                print(f"     Nie udało się wyrównać połączonych sekwencji")
        except Exception as e:
            print(f"     Błąd podczas alignowania połączonych sekwencji: {e}")
            import traceback
            traceback.print_exc()

        return residue_mapping, missing_atoms_report

    def create_residue_mapping(self, chain_mapping, target_residues, model_residues):
        """Tworzy mapowanie reszt na podstawie mapowania łańcuchów."""
        residue_mapping = {}
        missing_atoms_report = []
        
        # Sprawdź czy mamy przypadek połączonego łańcucha
        if len(chain_mapping) > 0:
            model_chains_used = set(chain_mapping.values())
            if len(model_chains_used) == 1 and len(chain_mapping) > 1:
                # Model ma jeden łańcuch mapowany na wiele łańcuchów targetu
                model_chain = list(model_chains_used)[0]
                target_chains = list(chain_mapping.keys())
                
                # Utwórz sekwencje tylko dla zmapowanych łańcuchów targetu
                target_sequences_for_mapping = {}
                target_residues_for_mapping = {}
                for target_chain in target_chains:
                    if target_chain in target_residues:
                        target_sequences_for_mapping[target_chain] = ''.join([r.resname.strip() for r in target_residues[target_chain]])
                        target_residues_for_mapping[target_chain] = target_residues[target_chain]
                
                return self.create_residue_mapping_for_combined_chain(
                    target_sequences_for_mapping, target_residues_for_mapping, 
                    model_chain, model_residues
                )
        
        # Normalne mapowanie łańcuchów
        for target_chain, model_chain in chain_mapping.items():
            if model_chain not in model_residues or target_chain not in target_residues:
                print(f"  OSTRZEŻENIE: Brak danych dla {model_chain}->{target_chain}")
                continue
            
            model_chain_residues = model_residues[model_chain]
            target_chain_residues = target_residues[target_chain]
            
            # Sprawdź czy sekwencje się zgadzają
            model_seq = ''.join([r.resname.strip() for r in model_chain_residues])
            target_seq = ''.join([r.resname.strip() for r in target_chain_residues])
            
            print(f"  Mapowanie łańcucha {model_chain}->{target_chain}:")
            print(f"    Target: {target_seq} (długość: {len(target_seq)}, reszty: {[r.id[1] for r in target_chain_residues]})")
            print(f"    Model:  {model_seq} (długość: {len(model_seq)}, reszty: {[r.id[1] for r in model_chain_residues]})")
            
            if model_seq == target_seq:
                # Sekwencje identyczne - mapuj jedna do jednej
                min_len = min(len(model_chain_residues), len(target_chain_residues))
                for i in range(min_len):
                    model_res = model_chain_residues[i]
                    target_res = target_chain_residues[i]
                    
                    # Sprawdź atomy
                    missing_atoms = self.check_missing_atoms(target_res, model_res)
                    if missing_atoms:
                        missing_atoms_report.append({
                            'model_chain': model_chain,
                            'model_residue': model_res.id[1],
                            'target_chain': target_chain,
                            'target_residue': target_res.id[1],
                            'missing_atoms': missing_atoms
                        })
                    
                    # Dodaj mapowanie
                    residue_mapping[(model_chain, model_res.id[1])] = (target_chain, target_res.id[1])
                    print(f"    {model_chain}{model_res.id[1]} -> {target_chain}{target_res.id[1]}")
                
                if len(model_chain_residues) > min_len:
                    print(f"    ! Pominięto {len(model_chain_residues) - min_len} reszt w modelu (za długie)")
                if len(target_chain_residues) > min_len:
                    print(f"    ! Pominięto {len(target_chain_residues) - min_len} reszt w targetcie (za długie)")
                    
            else:
                # Sekwencje różne - użyj alignowania
                print(f"    Sekwencje różne - używam alignowania")
                try:
                    alignments = self.aligner.align(target_seq, model_seq)
                    if alignments:
                        best_alignment = alignments[0]
                        aligned_target, aligned_model = best_alignment[0], best_alignment[1]
                        
                        print(f"    Wyrównanie target: {aligned_target}")
                        print(f"    Wyrównanie model:  {aligned_model}")
                        
                        # Mapowanie po alignmencie
                        target_idx, model_idx = 0, 0
                        mapped_count = 0
                        for t_char, m_char in zip(aligned_target, aligned_model):
                            if t_char != '-' and m_char != '-':
                                # Mapuj te reszty
                                if (target_idx < len(target_chain_residues) and 
                                    model_idx < len(model_chain_residues)):
                                    target_res = target_chain_residues[target_idx]
                                    model_res = model_chain_residues[model_idx]
                                    
                                    # Sprawdź czy nukleotydy się zgadzają
                                    if target_res.resname.strip() == model_res.resname.strip():
                                        # Sprawdź atomy
                                        missing_atoms = self.check_missing_atoms(target_res, model_res)
                                        if missing_atoms:
                                            missing_atoms_report.append({
                                                'model_chain': model_chain,
                                                'model_residue': model_res.id[1],
                                                'target_chain': target_chain,
                                                'target_residue': target_res.id[1],
                                                'missing_atoms': missing_atoms
                                            })
                                        
                                        # Dodaj mapowanie
                                        residue_mapping[(model_chain, model_res.id[1])] = (target_chain, target_res.id[1])
                                        print(f"    {model_chain}{model_res.id[1]} -> {target_chain}{target_res.id[1]} (po alignmencie)")
                                        mapped_count += 1
                                    else:
                                        print(f"     Różne nukleotydy: {model_chain}{model_res.id[1]}({model_res.resname}) != {target_chain}{target_res.id[1]}({target_res.resname})")
                                else:
                                    print(f"    ! Indeks poza zakresem: target_idx={target_idx}, model_idx={model_idx}")
                                
                                target_idx += 1
                                model_idx += 1
                            elif t_char == '-':
                                # Gap w targetcie - pomiń resztę w modelu
                                if model_idx < len(model_chain_residues):
                                    print(f"    - Gap w targetcie, pomijam {model_chain}{model_chain_residues[model_idx].id[1]}")
                                    model_idx += 1
                            elif m_char == '-':
                                # Gap w modelu - pomiń resztę w targetcie
                                if target_idx < len(target_chain_residues):
                                    print(f"    - Gap w modelu, pomijam {target_chain}{target_chain_residues[target_idx].id[1]}")
                                    target_idx += 1
                        
                        print(f"    Zmapowano {mapped_count} reszt po alignmencie")
                    else:
                        print(f"     Brak alignowania dla {target_chain} i {model_chain}")
                except Exception as e:
                    print(f"     Błąd podczas alignowania: {e}")
                    import traceback
                    traceback.print_exc()
        
        return residue_mapping, missing_atoms_report
    
    def check_missing_atoms(self, target_residue, model_residue):
        """Sprawdza czy w modelu brakuje atomów które są w targetcie."""
        missing_atoms = []
        
        # Lista atomów które powinny być w nukleotydzie RNA
        expected_atoms = ['P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'', 'O3\'', 'C2\'', 'O2\'', 'C1\'']
        
        # Dodaj atomy specyficzne dla bazy
        base = target_residue.resname.strip()
        if base == 'A':
            expected_atoms.extend(['N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4'])
        elif base == 'G':
            expected_atoms.extend(['N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4'])
        elif base == 'C':
            expected_atoms.extend(['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'])
        elif base == 'U':
            expected_atoms.extend(['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'])
        
        for atom_name in expected_atoms:
            try:
                target_residue[atom_name]  # Sprawdź czy atom istnieje w targetcie
                try:
                    model_residue[atom_name]  # Sprawdź czy atom istnieje w modelu
                except KeyError:
                    missing_atoms.append(atom_name)
            except KeyError:
                # Ten atom nie jest obecny w targetcie, więc nie ma co sprawdzać
                pass
        
        return missing_atoms
    
    def calculate_rmsd_with_mapping(self, target_structure, model_structure, mapping):
        """Oblicza RMSD używając mapowania."""
        if len(mapping) < 3:
            print("Za mało zmapowanych reszt do obliczenia RMSD (min. 3)")
            return float('inf')
        
        target_atoms = []
        model_atoms = []
        
        for model_full_id, target_full_id in mapping.items():
            model_chain_id, model_res_id = model_full_id
            target_chain_id, target_res_id = target_full_id
            
            # Znajdź atomy P w obu strukturach
            model_atom = self.find_atom_object(model_structure, model_chain_id, model_res_id, 'P')
            target_atom = self.find_atom_object(target_structure, target_chain_id, target_res_id, 'P')
            
            if model_atom is None or target_atom is None:
                # Spróbuj z atomem C4'
                model_atom = self.find_atom_object(model_structure, model_chain_id, model_res_id, 'C4\'')
                target_atom = self.find_atom_object(target_structure, target_chain_id, target_res_id, 'C4\'')
            
            if model_atom is not None and target_atom is not None:
                target_atoms.append(target_atom)
                model_atoms.append(model_atom)
        
        if len(target_atoms) < 3:
            print("Za mało atomów do obliczenia RMSD (min. 3)")
            return float('inf')
        
        # Oblicz RMSD
        from Bio.PDB import Superimposer
        superimposer = Superimposer()
        superimposer.set_atoms(target_atoms, model_atoms)
        rmsd = superimposer.rms
        
        return rmsd
    
    def find_atom_object(self, structure, chain_id, res_id, atom_name):
        """Znajduje obiekt Atom w strukturze."""
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    for residue in chain:
                        if residue.id[1] == res_id:
                            try:
                                return residue[atom_name]
                            except KeyError:
                                return None
        return None



    def fix_pdb_file(self, input_pdb, output_pdb, mapping, target_chain_order, remove_extra_atoms=False):
        """Naprawia plik PDB modelu, zmieniając łańcuchy, numery reszt i numerację atomów."""
        print(f"Naprawianie pliku PDB: {input_pdb} -> {output_pdb}")
        if remove_extra_atoms:
            print("  Tryb: usuwanie nadmiarowych atomów i reszt WŁĄCZONY")
        else:
            print("  Tryb: ZACHOWYWANIE wszystkich atomów w zmapowanych resztach")
        
        # Odczytaj oryginalny plik PDB
        with open(input_pdb, 'r') as f:
            lines = f.readlines()
        
        # Znajdź indeksy first_atom_index i last_atom_index
        first_atom_index = None
        last_atom_index = None
        for i, line in enumerate(lines):
            if line.startswith(('ATOM', 'HETATM')):
                if first_atom_index is None:
                    first_atom_index = i
                last_atom_index = i
        
        # Jeśli nie ma atomów, to zwróć oryginał
        if first_atom_index is None:
            with open(output_pdb, 'w') as f:
                f.writelines(lines)
            print("Brak atomów w pliku PDB.")
            return 0, 0, 0
        
        # Linie przed atomami
        before_atoms = lines[:first_atom_index]
        # Linie po atomach
        after_atoms = lines[last_atom_index+1:]
        
        # Przetwarzanie linii atomów
        processed_atom_lines = []
        
        unmapped_count = 0
        mapped_count = 0
        removed_atoms_count = 0
        removed_residues_count = 0
        
        current_residue = None
        current_residue_has_mapping = False
        
        for i in range(first_atom_index, last_atom_index+1):
            line = lines[i]
            if line.startswith(('ATOM', 'HETATM')):
                chain_id = line[21]
                try:
                    resid = int(line[22:26].strip())
                except ValueError:
                    # Jeśli nie można sparsować numeru reszty, pomiń
                    unmapped_count += 1
                    continue
                
                # Sprawdź czy zmieniła się reszta
                if current_residue != (chain_id, resid):
                    current_residue = (chain_id, resid)
                    current_residue_has_mapping = (chain_id, resid) in mapping
                    if not current_residue_has_mapping:
                        removed_residues_count += 1
                
                if current_residue_has_mapping:
                    # Reszta jest zmapowana - ZACHOWUJEMY WSZYSTKIE ATOMY
                    new_chain, new_resid = mapping[(chain_id, resid)]
                    
                    # Bez usuwania atomów - dodaj wszystkie atomy zmapowanych reszt
                    new_line = line[:21] + new_chain + line[22:]
                    new_resid_str = str(new_resid).rjust(4)
                    new_line = new_line[:22] + new_resid_str + new_line[26:]
                    processed_atom_lines.append(new_line)
                    mapped_count += 1
                else:
                    # Reszta nie jest zmapowana - pomiń
                    unmapped_count += 1
        
        print(f"  Zmapowane atomy: {mapped_count}")
        print(f"  Usunięte reszty: {removed_residues_count}")
        print(f"  Niezmapowane atomy: {unmapped_count}")
        
        # Posortuj atomy: najpierw według kolejności łańcuchów w targetcie, potem według numeru reszty
        atom_lines_dict = defaultdict(list)
        
        for line in processed_atom_lines:
            chain_id = line[21]
            try:
                resid = int(line[22:26].strip())
            except ValueError:
                continue
            atom_lines_dict[(chain_id, resid)].append(line)
        
        sorted_atom_lines = []
        
        # Łańcuchy w kolejności targetu
        for chain_id in target_chain_order:
            # Pobierz wszystkie reszty dla tego łańcucha
            residues = [ (c, r) for (c, r) in atom_lines_dict.keys() if c == chain_id ]
            residues.sort(key=lambda x: x[1])
            for (c, r) in residues:
                sorted_atom_lines.extend(atom_lines_dict[(c, r)])
        
        # Łańcuchy nieobecne w targetcie (powinny być rzadkie)
        remaining_chains = set([c for (c, r) in atom_lines_dict.keys()]) - set(target_chain_order)
        for chain_id in sorted(remaining_chains):  # sortujemy alfabetycznie
            residues = [ (c, r) for (c, r) in atom_lines_dict.keys() if c == chain_id ]
            residues.sort(key=lambda x: x[1])
            for (c, r) in residues:
                sorted_atom_lines.extend(atom_lines_dict[(c, r)])
        
        # NUMERUJ ATOMY OD 1 - NAPRAWIONA NUMERACJA
        atom_num = 1
        new_sorted_atom_lines = []
        for line in sorted_atom_lines:
            # Zamień numer atomu (kolumny 7-11) na nowy numer
            # Format PDB: kolumny 7-11 to pozycje 6-10 (0-based indexing)
            atom_num_str = str(atom_num).rjust(5)
            new_line = line[:6] + atom_num_str + line[11:]
            new_sorted_atom_lines.append(new_line)
            atom_num += 1
        
        # Połącz wszystko
        new_lines = before_atoms + new_sorted_atom_lines + after_atoms
        
        # Zapisz
        with open(output_pdb, 'w') as f:
            f.writelines(new_lines)
        
        changed_count = len(new_sorted_atom_lines)
        print(f"Naprawiono {changed_count} atomów w pliku PDB (łańcuchy, reszty i numery atomów).")
        return changed_count, 0, removed_residues_count  # removed_atoms_count zawsze 0

    def create_trimmed_solution_pdb(self, target_pdb, output_solution_pdb, residue_mapping, model_structure=None):
        """Tworzy przycięty plik solution zawierający tylko zmapowane reszty i atomy."""
        print(f"Tworzenie przyciętego solution: {target_pdb} -> {output_solution_pdb}")
        
        # Wczytaj strukturę targetu
        parser = PDBParser()
        target_structure = parser.get_structure("target", target_pdb)
        
        # Odwróć mapowanie: (target_chain, target_resid) -> (model_chain, model_resid)
        reverse_mapping = {v: k for k, v in residue_mapping.items()}
        
        # Zbierz wszystkie zmapowane reszty targetu
        target_residues_to_keep = set(reverse_mapping.keys())
        
        # Odczytaj oryginalny plik target PDB
        with open(target_pdb, 'r') as f:
            lines = f.readlines()
        
        # Znajdź indeksy first_atom_index i last_atom_index
        first_atom_index = None
        last_atom_index = None
        for i, line in enumerate(lines):
            if line.startswith(('ATOM', 'HETATM')):
                if first_atom_index is None:
                    first_atom_index = i
                last_atom_index = i
        
        # Jeśli nie ma atomów, to zwróć oryginał
        if first_atom_index is None:
            with open(output_solution_pdb, 'w') as f:
                f.writelines(lines)
            print("Brak atomów w pliku target PDB.")
            return 0
        
        # Linie przed atomami
        before_atoms = lines[:first_atom_index]
        # Linie po atomach
        after_atoms = lines[last_atom_index+1:]
        
        # Przetwarzanie linii atomów - zachowujemy tylko atomy które istnieją w modelu
        processed_atom_lines = []
        kept_residues = set()
        
        current_residue = None
        current_residue_should_keep = False
        current_model_residue = None
        
        for i in range(first_atom_index, last_atom_index+1):
            line = lines[i]
            if line.startswith(('ATOM', 'HETATM')):
                chain_id = line[21]
                try:
                    resid = int(line[22:26].strip())
                except ValueError:
                    continue
                
                # Sprawdź czy zmieniła się reszta
                if current_residue != (chain_id, resid):
                    current_residue = (chain_id, resid)
                    current_residue_should_keep = (chain_id, resid) in target_residues_to_keep
                    
                    # Znajdź odpowiadającą resztę w modelu
                    current_model_residue = None
                    if current_residue_should_keep and model_structure:
                        kept_residues.add(current_residue)
                        # Znajdź modelową resztę poprzez reverse mapping
                        if (chain_id, resid) in reverse_mapping:
                            model_chain_id, model_res_id = reverse_mapping[(chain_id, resid)]
                            # Znajdź resztę w modelu
                            for model in model_structure:
                                if model_chain_id in model:
                                    for res in model[model_chain_id]:
                                        if res.id[1] == model_res_id:
                                            current_model_residue = res
                                            break
                                if current_model_residue:
                                    break
                
                if current_residue_should_keep:
                    if model_structure and current_model_residue:
                        # Sprawdź czy ten atom istnieje w modelu
                        atom_name = line[12:16].strip()
                        
                        # Sprawdź czy atom istnieje w modelowej reszcie
                        atom_exists_in_model = False
                        try:
                            if current_model_residue[atom_name]:
                                atom_exists_in_model = True
                        except KeyError:
                            atom_exists_in_model = False
                        
                        # Tylko zachowaj atomy które istnieją w modelu
                        if atom_exists_in_model:
                            processed_atom_lines.append(line)
                    else:
                        # Jeśli nie mamy model_structure, zachowujemy wszystkie atomy
                        processed_atom_lines.append(line)
        
        print(f"  Zachowane reszty w solution: {len(kept_residues)}")
        print(f"  Zachowane atomy w solution: {len(processed_atom_lines)}")
        
        # NUMERUJ ATOMY OD 1
        atom_num = 1
        new_sorted_atom_lines = []
        for line in processed_atom_lines:
            atom_num_str = str(atom_num).rjust(5)
            new_line = line[:6] + atom_num_str + line[11:]
            new_sorted_atom_lines.append(new_line)
            atom_num += 1
        
        # Połącz wszystko
        new_lines = before_atoms + new_sorted_atom_lines + after_atoms
        
        # Zapisz
        with open(output_solution_pdb, 'w') as f:
            f.writelines(new_lines)
        
        changed_count = len(new_sorted_atom_lines)
        print(f"Utworzono przycięty solution z {changed_count} atomami.")
        return changed_count

    def verify_fixed_pdb(self, target_pdb, fixed_pdb, residue_mapping):
        """Sprawdza atomy w naprawionym pliku PDB względem targetu."""
        print(f"\nWeryfikacja atomów w naprawionym pliku PDB...")
        
        # Wczytaj struktury
        parser = PDBParser()
        target_structure = parser.get_structure("target", target_pdb)
        fixed_structure = parser.get_structure("fixed", fixed_pdb)
        
        # Odwróć mapowanie: (target_chain, target_resid) -> (model_chain, model_resid)
        reverse_mapping = {v: k for k, v in residue_mapping.items()}
        
        atom_comparison_log = []
        total_extra_atoms = 0
        total_missing_atoms = 0
        
        for target_chain_id, target_res_id in reverse_mapping.keys():
            # Znajdź resztę w targetcie
            target_residue = None
            for residue in target_structure[0][target_chain_id]:
                if residue.id[1] == target_res_id:
                    target_residue = residue
                    break
            
            # Znajdź odpowiadającą resztę w naprawionym modelu
            model_chain_id, model_res_id = reverse_mapping[(target_chain_id, target_res_id)]
            fixed_residue = None
            for residue in fixed_structure[0][target_chain_id]:
                if residue.id[1] == model_res_id:
                    fixed_residue = residue
                    break
            
            if target_residue and fixed_residue:
                # Zbierz atomy z targetu
                target_atoms = set()
                for atom in target_residue:
                    target_atoms.add(atom.name)
                
                # Zbierz atomy z naprawionego modelu
                fixed_atoms = set()
                for atom in fixed_residue:
                    fixed_atoms.add(atom.name)
                
                # Znajdź różnice
                extra_atoms = fixed_atoms - target_atoms  # Atomy w modelu, których nie ma w targetcie
                missing_atoms = target_atoms - fixed_atoms  # Atomy w targetcie, których nie ma w modelu
                
                if extra_atoms or missing_atoms:
                    # Dodaj do logu
                    for atom_name in sorted(extra_atoms):
                        atom_comparison_log.append(f"+ ATOM {atom_name:>4} {fixed_residue.resname} {model_chain_id} {model_res_id:>3}")
                        total_extra_atoms += 1
                    
                    for atom_name in sorted(missing_atoms):
                        atom_comparison_log.append(f"- ATOM {atom_name:>4} {target_residue.resname} {target_chain_id} {target_res_id:>3}")
                        total_missing_atoms += 1
        
        return atom_comparison_log, total_extra_atoms, total_missing_atoms


    def process_all(self, target_file, model_file, output_pdb, remove_extra_atoms=False):
        """Główna metoda przetwarzająca wszystko."""
        print(f"=== DEBUG: Rozpoczynam process_all ===")
        print(f"  target_file: {target_file}")
        print(f"  model_file: {model_file}")
        print(f"  output_pdb: {output_pdb}")
        
        try:
            print("Ładowanie struktur...")
            parser = PDBParser()
            target_structure = parser.get_structure("target", target_file)
            model_structure = parser.get_structure("model", model_file)
            print(f"  DEBUG: Target structure loaded, models: {len(list(target_structure))}")
            print(f"  DEBUG: Model structure loaded, models: {len(list(model_structure))}")
            
            print("Ekstrakcja sekwencji i reszt...")
            target_sequences, target_residues = self.extract_chain_sequences(target_structure)
            model_sequences, model_residues = self.extract_chain_sequences(model_structure)
            
            print(f"  DEBUG: Target sequences: {target_sequences}")
            print(f"  DEBUG: Model sequences: {model_sequences}")
            print(f"  DEBUG: Target residues keys: {list(target_residues.keys())}")
            print(f"  DEBUG: Model residues keys: {list(model_residues.keys())}")
            
            # Ustaw atrybuty PRZED wywołaniem find_chain_mapping
            self._last_target_residues = target_residues
            self._last_model_residues = model_residues
            
            # Policz nukleotydy
            target_rna_count = self.count_rna_residues(target_structure)
            model_rna_count = self.count_rna_residues(model_structure)
            
            print(f"Target: {len(target_sequences)} łańcuchów, {target_rna_count} nukleotydów")
            print(f"Model:  {len(model_sequences)} łańcuchów, {model_rna_count} nukleotydów")
            
            print(f"\nMapowanie łańcuchów...")
            print(f"  DEBUG: Wywołuję find_chain_mapping")
            chain_mapping = self.find_chain_mapping(target_sequences, model_sequences)
            
            print(f"  DEBUG: chain_mapping = {chain_mapping}")
            
            if not chain_mapping:
                print("Nie znaleziono pasujących łańcuchów!")
                return {}, [], float('inf'), 0, 0, 0, [], 0, 0, None
            
            print(f"\nMapowanie reszt i sprawdzanie atomów...")
            residue_mapping, missing_atoms_report = self.create_residue_mapping(
                chain_mapping, target_residues, model_residues
            )
            
            print(f"Zmapowano {len(residue_mapping)} reszt")
            
            if missing_atoms_report:
                print(f"\nZNALEZIONO BRAKUJĄCE ATOMY ({len(missing_atoms_report)} reszt):")
                for report in missing_atoms_report:
                    print(f"  {report['model_chain']}{report['model_residue']} -> {report['target_chain']}{report['target_residue']}: brakuje {report['missing_atoms']}")
            else:
                print("\nWSZYSTKIE ATOMY ZGODNE!")
            
            # Oblicz RMSD
            rmsd = self.calculate_rmsd_with_mapping(target_structure, model_structure, residue_mapping)
            print(f"\nRMSD dla zmapowanych reszt: {rmsd:.3f} A")
            
            # Pobierz kolejność łańcuchów z targetu
            target_chain_order = self.get_target_chain_order(target_structure)
            print(f"Kolejność łańcuchów w targetcie: {target_chain_order}")
            
            # Napraw plik PDB modelu
            print(f"\nNaprawianie pliku PDB...")
            changed_count, removed_atoms_count, removed_residues_count = self.fix_pdb_file(
                model_file, output_pdb, residue_mapping, target_chain_order, remove_extra_atoms
            )
            
            print(f"  DEBUG: fix_pdb_file zakończone: changed={changed_count}, removed_residues={removed_residues_count}")
            
            # WERYFIKUJ naprawiony plik PDB
            atom_comparison_log, total_extra_atoms, total_missing_atoms = self.verify_fixed_pdb(
                target_file, output_pdb, residue_mapping
            )
            
            print(f"  DEBUG: verify_fixed_pdb: total_extra={total_extra_atoms}, total_missing={total_missing_atoms}")
            
            # Tworzenie pliku solution jeśli są brakujące atomy/reszty
            solution_output = None
            
            missing_residues_in_model = (len(residue_mapping) < target_rna_count)
            should_create_solution = missing_residues_in_model or total_missing_atoms > 0
            
            if should_create_solution:
                model_basename = os.path.basename(model_file)
                model_name = os.path.splitext(model_basename)[0]
                if model_name.startswith('PZ2_'):
                    model_suffix = model_name[4:]
                else:
                    model_suffix = model_name
                
                output_dir = os.path.dirname(output_pdb)
                solution_filename = f"solution_{model_suffix}.pdb"
                
                if output_dir:
                    solution_output = os.path.join(output_dir, solution_filename)
                else:
                    solution_output = solution_filename
                
                if output_dir and not os.path.exists(output_dir):
                    os.makedirs(output_dir, exist_ok=True)
                
                print(f"  DEBUG: Tworzę solution: {solution_output}")
                
                # Użyj poprawionej wersji create_trimmed_solution_pdb
                solution_atom_count = self.create_trimmed_solution_pdb(
                    target_file, solution_output, residue_mapping, model_structure
                )
                print(f"Utworzono plik solution z {solution_atom_count} atomami")
            else:
                print("Brak BRAKUJĄCYCH atomów/reszt w modelu - plik solution nie został utworzony")
            
            print(f"\n=== DEBUG: process_all zakończone ===")
            print(f"  residue_mapping: {len(residue_mapping)} wpisów")
            print(f"  solution_output: {solution_output}")
            
            return residue_mapping, missing_atoms_report, rmsd, changed_count, removed_atoms_count, removed_residues_count, atom_comparison_log, total_extra_atoms, total_missing_atoms, solution_output
        
        except Exception as e:
            print(f"\n=== BŁĄD W process_all ===")
            print(f"Typ błędu: {type(e).__name__}")
            print(f"Wiadomość: {e}")
            import traceback
            traceback.print_exc()
            return {}, [], float('inf'), 0, 0, 0, [], 0, 0, None
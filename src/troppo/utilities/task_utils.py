
import pandas as pd
import re
from cobra.core import Model

class TextTaskIO:
    """
    Parses RAVEN-style metabolic task text files (tab-separated).
    """
    def __init__(self, file_path):
        self.file_path = file_path
        self.tasks = self._parse_file()

    def _parse_file(self):
        tasks = []
        with open(self.file_path, 'r', encoding='utf-8') as f:
            lines = [l.rstrip() for l in f.readlines()]

        current_task = None
        global_offset = 0 
        offset_set = False

        for line in lines:
            if line.startswith("#"):
                continue
            if not line.strip():
                continue

            parts = line.split('\t')
            
            # Determine offset from the first task line we encounter
            if not offset_set:
                if parts[0].strip():
                    global_offset = 0
                    offset_set = True
                elif len(parts) > 1 and parts[1].strip():
                    global_offset = 1
                    offset_set = True
            
            # Identify columns based on offset
            # ID: 0 + offset
            # Desc: 1 + offset
            # ShouldFail: 2 + offset
            # IN: 3 + offset
            # OUT: 6 + offset
            
            idx_id = 0 + global_offset
            idx_desc = 1 + global_offset
            idx_sf = 2 + global_offset
            idx_in = 3 + global_offset
            idx_out = 6 + global_offset
            
            vid = parts[idx_id].strip() if len(parts) > idx_id else ""
            desc = parts[idx_desc].strip() if len(parts) > idx_desc else ""
            should_fail = parts[idx_sf].strip() if len(parts) > idx_sf else ""
            
            # IN and OUT might be in further columns
            inputs = parts[idx_in].strip() if len(parts) > idx_in else ""
            outputs = parts[idx_out].strip() if len(parts) > idx_out else ""
            
            if vid:
                # New Task
                current_task = {
                    'id': vid,
                    'description': desc,
                    'should_fail': str(should_fail).lower() == 'true', # Convert to bool carefully
                    'inputs': [],
                    'outputs': []
                }
                tasks.append(current_task)
            
            if current_task:
                # Add inputs
                if inputs:
                    # Split by ';'
                    mets = [x.strip() for x in inputs.split(';') if x.strip()]
                    current_task['inputs'].extend(mets)
                
                # Add outputs (optional, for objective)
                if outputs:
                    mets = [x.strip() for x in outputs.split(';') if x.strip()]
                    current_task['outputs'].extend(mets)

        return tasks

    def get_task_by_name(self, name_snippet):
        for t in self.tasks:
            if name_snippet.lower() in t['description'].lower():
                return t
        return None

def map_metabolites_to_exchanges(model: Model, metabolite_names: list, verbose=True):
    """
    Maps a list of metabolite names (e.g. 'arginine[e]') to Exchange Reactions in the model.
    Returns:
        dict: {metabolite_name: exchange_reaction_id}
        list: missing_metabolites
    """
    mapping = {}
    missing = []
    
    # Pre-cache model metabolites for speed
    # Map (name, compartment) -> id
    # Also (id) -> id
    
    met_map = {}
    for m in model.metabolites:
        # Key by clean ID
        met_map[m.id] = m
        # Key by name (lowercase) + compartment
        key = (m.name.lower(), m.compartment)
        if key not in met_map:
            met_map[key] = []
        met_map[key].append(m)
        
    for met_str in metabolite_names:
        # Parse "name[comp]" or just "name"
        # Regex for name[comp]
        match = re.match(r"(.+)\[([a-z])\]$", met_str)
        if match:
            name, comp = match.groups()
        else:
            name = met_str
            comp = 'e' # Default to extracellular if not specified?
            
        name_clean = name.strip()
        
        # Strategies to find metabolite:
        found_met = None
        
        # 1. Exact ID match (check if met_str is an ID)
        if met_str in model.metabolites:
            found_met = model.metabolites.get_by_id(met_str)
        
        # 2. Name + Compartment match
        if not found_met:
            # Prepare candidate keys
            candidate_keys = []
            candidate_keys.append((name_clean.lower(), comp))
            
            # Helper for HumanGEM: [e] usually maps to 's' compartment
            if comp == 'e':
                candidate_keys.append((name_clean.lower(), 's'))
            
            if verbose and "arginine" in name_clean:
                print(f"DEBUG: Tracing '{met_str}'...")
                print(f"  Regex parsed: name='{name_clean}', comp='{comp}'")
                print(f"  Candidate keys: {candidate_keys}")
                for k in candidate_keys:
                    print(f"  Checking key {k}: {'FOUND' if k in met_map else 'NOT FOUND'}")
            
            for key in candidate_keys:
                if key in met_map:
                    candidates = met_map[key]
                    if candidates:
                        found_met = candidates[0]
                        break
        
        # 3. Try generic naming (L-arginine vs arginine)
        if not found_met and not name_clean.lower().startswith('l-'):
            l_name = f"l-{name_clean.lower()}"
            candidate_keys = [(l_name, comp)]
            if comp == 'e':
                candidate_keys.append((l_name, 's'))
                
            for key in candidate_keys:
                if key in met_map:
                    found_met = met_map[key][0]
                    break
                    
        # Find Exchange Reaction
        if found_met:
            ex_rxn = None
            # Check reactions of metabolite
            for rxn in found_met.reactions:
                # Is it an exchange?
                # 1 product, 0 reactants OR 1 reactant, 0 products (boundary)
                if len(rxn.metabolites) == 1:
                    ex_rxn = rxn
                    break
                # Or revesible exchange A <=> (boundary)
            
            if ex_rxn:
                mapping[met_str] = ex_rxn.id
            else:
                if verbose: print(f"Warning: Metabolite '{found_met.id}' found but no exchange reaction.")
                missing.append(met_str)
        else:
            if verbose: print(f"Warning: Metabolite '{met_str}' not found in model.")
            missing.append(met_str)
            
    return mapping, missing

def get_task_medium(model: Model, task: dict):
    """
    Extracts medium composition (exchange reactions) from a task's inputs.
    Returns:
        list: List of exchange reaction IDs to be set as medium (open uptake).
        list: Missing inputs (metabolites not mapped).
    """
    mapping, missing = map_metabolites_to_exchanges(model, task['inputs'], verbose=False)
    return list(mapping.values()), missing

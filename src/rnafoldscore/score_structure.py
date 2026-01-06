import argparse

def main():
    p = argparse.ArgumentParser(
        prog="rnafoldscore.score_structure",
        description="Score a PDB structure using trained profiles (skeleton)."
    )
    p.add_argument("--pdb", required=True)
    p.add_argument("--profiles", required=True)
    args = p.parse_args()

    print("SCORE skeleton OK")
    print("pdb:", args.pdb)
    print("profiles:", args.profiles)

if __name__ == "__main__":
    main()
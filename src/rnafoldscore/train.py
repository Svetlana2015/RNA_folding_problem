import argparse

def main():
    p = argparse.ArgumentParser(
        prog="rnafoldscore.train",
        description="Train distance-based profiles from PDB directory (skeleton)."
    )
    p.add_argument("--pdb-dir", required=True)
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()

    print("TRAIN skeleton OK")
    print("pdb_dir:", args.pdb_dir)
    print("out_dir:", args.out_dir)

if __name__ == "__main__":
    main()
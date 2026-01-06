import argparse

def main():
    p = argparse.ArgumentParser(
        prog="rnafoldscore.plot_profiles",
        description="Plot profiles (skeleton)."
    )
    p.add_argument("--profiles", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    print("PLOT skeleton OK")
    print("profiles:", args.profiles)
    print("out:", args.out)

if __name__ == "__main__":
    main()
    
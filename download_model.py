"""Download the e_coli_core metabolic model using COBRApy."""

import cobra
from pathlib import Path


def main():
    data_dir = Path(__file__).parent / "data"
    data_dir.mkdir(exist_ok=True)

    output_path = data_dir / "e_coli_core.xml"

    if output_path.exists():
        print(f"Model already exists at {output_path}")
    else:
        print("Downloading e_coli_core model...")
        model = cobra.io.load_model("textbook")  # e_coli_core is the "textbook" model
        cobra.io.write_sbml_model(model, str(output_path))
        print(f"Model saved to {output_path}")

    # Verify by loading it back
    model = cobra.io.read_sbml_model(str(output_path))
    print(f"Verification: loaded model '{model.id}' with "
          f"{len(model.reactions)} reactions, "
          f"{len(model.metabolites)} metabolites, "
          f"{len(model.genes)} genes.")


if __name__ == "__main__":
    main()

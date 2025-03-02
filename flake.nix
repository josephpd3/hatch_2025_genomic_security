{
  description = "Python development environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };

        # Specify the Python version you want to use
        python = pkgs.python311;

        # Create a Python environment with the packages you need
        pythonEnv = python.withPackages (ps: with ps; [
          # Add the Python packages you need
          numpy
          polars
          requests
          black
          biopython

          # GRPC deps
          grpcio-tools
        ]);
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = with pkgs; [
            # Your custom Python environment
            pythonEnv

            go
            gopls
            go-tools
            delve

            gotest

            # Additional development tools
            pkgs.git
            pkgs.pre-commit
          ];

          # Environment variables
          shellHook = ''
            # Create a .venv directory that points to our nix Python environment
            # This helps tools like VSCode find the right Python environment
            if [ ! -e .venv ]; then
              ln -s ${pythonEnv} .venv
            fi

            # Tells pip to put packages into $PIP_PREFIX instead of the nix store
            export PIP_PREFIX=$(pwd)/_build/pip_packages
            export PYTHONPATH="$PIP_PREFIX/${python.sitePackages}:$PYTHONPATH"
            export PATH="$PIP_PREFIX/bin:$PATH"

            # Setup pre-commit if a config exists
            if [ -f .pre-commit-config.yaml ]; then
              pre-commit install
            fi

            echo "Python development environment activated!"
          '';
        };
      }
    );
}
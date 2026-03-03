###############################################################################
#                         Multi-Package Source Exporter                       #
###############################################################################

# Dynamically identify all subdirectories (each is treated as a package)
PACKAGES := $(patsubst %/,%,$(wildcard */))

# File extensions to collect (Mathematica script and package files)
EXTENSIONS = -name "*.m" -o -name "*.wl"

# Default target: export sources for all detected packages
.PHONY: all
all: export-all

# Target to export all packages in one command
.PHONY: export-all $(PACKAGES)
export-all: $(PACKAGES)

# ---------------------------------------------------------------------------
# Pattern rule to export a specific package.
# Usage: make <package_name> (e.g., make QuantumWalks)
#
# Concatenates all .m and .wl files into a single text file.
# Each file is preceded by a header showing its path within the package.
# Files are sorted to ensure deterministic output.
# ---------------------------------------------------------------------------
$(PACKAGES):
	@echo "Exporting source for package: $@"
	@find $@ -type f \( $(EXTENSIONS) \) | sort | while read -r file; do \
		echo "=== FILENAME: $${file} ==="; \
		cat "$$file"; \
		printf "\n\n"; \
	done > $@_Source.txt
	@echo "$@_Source.txt generated successfully."


###############################################################################
#                                 CLEAN SECTION                               #
###############################################################################

# ---------------------------------------------------------------------------
# clean:
#   Removes all generated *_Source.txt files from the root directory.
# ---------------------------------------------------------------------------

.PHONY: clean

clean:
	@echo "Cleaning up exported source files..."
	rm -f *_Source.txt
	@echo "Cleanup complete."

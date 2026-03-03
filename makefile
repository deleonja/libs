###############################################################################
#                          MULTI-PACKAGE EXPORT SECTION                       #
###############################################################################

# Dynamically identify all subdirectories (each is treated as a package)
PACKAGES := $(patsubst %/,%,$(wildcard */))

# File extensions to collect (Mathematica script and package files)
EXTENSIONS = -name "*.m" -o -name "*.wl"

.PHONY: all export-all $(PACKAGES)

# Default target: export sources for all detected packages
all: export-all

# Target to export all packages in one command
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
#                              VERSIONING SECTION                             #
###############################################################################

.PHONY: bump

# ---------------------------------------------------------------------------
# bump:
#   Increments the semantic version (Major.Minor.Patch) in version.txt
#   and updates the Version field in PacletInfo.wl.
#   Automatically stages and commits these changes to the current branch.
#
#   Usage:
#     make bump type=major
#     make bump type=minor
#     make bump type=patch
# ---------------------------------------------------------------------------
bump:
ifndef type
	$(error You must specify a type. Example: make bump type=patch)
endif
	@CurrentVersion=$$(cat version.txt); \
	Major=$$(echo $$CurrentVersion | cut -d. -f1); \
	Minor=$$(echo $$CurrentVersion | cut -d. -f2); \
	Patch=$$(echo $$CurrentVersion | cut -d. -f3); \
	if [ "$(type)" = "major" ]; then \
		Major=$$((Major + 1)); Minor=0; Patch=0; \
	elif [ "$(type)" = "minor" ]; then \
		Minor=$$((Minor + 1)); Patch=0; \
	elif [ "$(type)" = "patch" ]; then \
		Patch=$$((Patch + 1)); \
	else \
		echo "Error: Define type as 'major', 'minor', or 'patch'"; \
		exit 1; \
	fi; \
	NewVersion="$$Major.$$Minor.$$Patch"; \
	echo $$NewVersion > version.txt; \
	sed -i "s/^\([ \t]*\)Version -> \"[^\"]*\"/\1Version -> \"$$NewVersion\"/" PacletInfo.wl; \
	git add version.txt PacletInfo.wl; \
	git commit -m "Bump version to $$NewVersion"; \
	echo "Successfully bumped $(type) version and committed. New version: $$NewVersion"


###############################################################################
#                              GIT RELEASE SECTION                            #
###############################################################################

# Dynamically identify the package name from the current directory
PKG_NAME ?= $(shell basename $(CURDIR))

.PHONY: release

# ---------------------------------------------------------------------------
# release:
#   Automates the git release workflow:
#   1. Validates current directory (QMB or QuantumWalks).
#   2. Merges develop into main and creates a version tag.
#   3. Pushes changes to origin with retry logic for SSH stability.
#   4. Merges main back into develop to keep branches synchronized.
# ---------------------------------------------------------------------------
release:
	@if [ "$(PKG_NAME)" != "QMB" ] && [ "$(PKG_NAME)" != "QuantumWalks" ]; then \
		echo "Error: Release can only be run from the QMB or QuantumWalks directories."; \
		exit 1; \
	fi; \
	CurrentVersion=$$(cat version.txt); \
	echo "Starting release process for $(PKG_NAME) v$$CurrentVersion..."; \
	git checkout main || { echo "Error: Failed to checkout main."; exit 1; }; \
	git merge develop -m "Merge develop into main for release $(PKG_NAME)-v$$CurrentVersion" || { echo "Error: Merge failed."; exit 1; }; \
	git tag -a "$(PKG_NAME)-v$$CurrentVersion" -m "Release $(PKG_NAME)-v$$CurrentVersion"; \
	echo "Pushing to remote..."; \
	sleep 1; \
	git push origin main --tags || { echo "Retrying push..."; sleep 2; git push origin main --tags; } || { echo "Error: Push failed."; exit 1; }; \
	echo "Synchronizing develop with main..."; \
	git checkout develop; \
	git merge main -m "Sync develop with main after release $(PKG_NAME)-v$$CurrentVersion"; \
	git push origin develop; \
	echo "Release $(PKG_NAME)-v$$CurrentVersion successfully merged, tagged, and synchronized."


###############################################################################
#                                CLEAN SECTION                                #
###############################################################################

.PHONY: clean

# ---------------------------------------------------------------------------
# clean:
#   Removes all generated *_Source.txt files from the root directory.
# ---------------------------------------------------------------------------
clean:
	@echo "Cleaning up exported source files..."
	@rm -f *_Source.txt
	@echo "Cleanup complete."

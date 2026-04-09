###############################################################################
#                          MULTI-PACKAGE EXPORT SECTION                       #
###############################################################################

# Dynamically identify all subdirectories (each is treated as a package)
PACKAGES := $(patsubst %/,%,$(wildcard */))

# File extensions to collect (Mathematica script and package files)
EXTENSIONS = -name "*.m" -o -name "*.wl"

# Current git branch
BRANCH := $(shell git rev-parse --abbrev-ref HEAD)

.PHONY: all export-all $(PACKAGES)

# Default target: export sources for all detected packages
all: export-all

# Target to export all packages in one command
export-all: $(PACKAGES)

$(PACKAGES):
	@echo "Exporting source for package: $@"
	@VER=$$(cat $@/version.txt); \
	find $@ -type f \( $(EXTENSIONS) \) | sort | while read -r file; do \
		echo "=== FILENAME: $${file} ==="; \
		cat "$$file"; \
		printf "\n\n"; \
	done > "$@_v$${VER}_Source.txt"; \
	echo "$@_v$${VER}_Source.txt generado con éxito."


###############################################################################
#                               VERSIONING SECTION                            #
###############################################################################

# Detectar el Sistema Operativo para compatibilidad de 'sed'
UNAME_S := $(shell uname -s)
# Obtener el nombre del paquete actual (ej. QMB o QuantumWalks)
CURRENT_PKG := $(shell basename $(CURDIR))

.PHONY: bump

bump:
ifndef type
	$(error You must specify a type. Example: make bump type=patch)
endif
	@if [ "$(BRANCH)" != "develop" ]; then \
		echo "Error: Version bumping must be performed on the 'develop' branch."; \
		exit 1; \
	fi
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
	echo "Updating files to v$$NewVersion..."; \
	if [ "$(UNAME_S)" = "Darwin" ]; then \
		sed -i '' "s/^\([ \t]*\)Version -> \"[^\"]*\"/\1Version -> \"$$NewVersion\"/" PacletInfo.wl; \
		sed -i '' "s/\(\*\*Current Version\*\*:\) v[0-9.]*/\1 v$$NewVersion/" README.md; \
		sed -i '' "s/\(\*\*\*$(CURRENT_PKG)\*\*\*:\) v[0-9.]*/\1 v$$NewVersion/" ../README.md; \
	else \
		sed -i "s/^\([ \t]*\)Version -> \"[^\"]*\"/\1Version -> \"$$NewVersion\"/" PacletInfo.wl; \
		sed -i "s/\(\*\*Current Version\*\*:\) v[0-9.]*/\1 v$$NewVersion/" README.md; \
		sed -i "s/\(\*\*\*$(CURRENT_PKG)\*\*\*:\) v[0-9.]*/\1 v$$NewVersion/" ../README.md; \
	fi; \
	git add version.txt PacletInfo.wl README.md ../README.md; \
	git commit -m "Bump $(CURRENT_PKG) version to $$NewVersion"; \
	echo "Successfully bumped $(type) version. New version: $$NewVersion"; \
	$(MAKE) release


###############################################################################
#                               GIT RELEASE SECTION                           #
###############################################################################

PKG_NAME ?= $(shell basename $(CURDIR))

.PHONY: release

release:
	@if [ "$(PKG_NAME)" != "QMB" ] && [ "$(PKG_NAME)" != "QuantumWalks" ]; then \
		echo "Error: Release can only be run from the QMB or QuantumWalks directories."; \
		exit 1; \
	fi; \
	if [ "$(BRANCH)" != "develop" ]; then \
		echo "Error: The release process must start from the 'develop' branch."; \
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
#                                  CLEAN SECTION                              #
###############################################################################

.PHONY: clean

clean:
	@echo "Cleaning up exported source files..."
	@rm -f *_v*_Source.txt
	@echo "Cleanup complete."

###############################################################################
#                          GOOGLE DRIVE UPLOAD SECTION                        #
###############################################################################

DriveRemote := gdrive
DrivePath := WolframSources

.PHONY: upload

# ---------------------------------------------------------------------------
# upload:
#   Diseñado para ejecutarse DESDE LA RAÍZ del proyecto.
# ---------------------------------------------------------------------------
upload:
ifndef pkg
	$(error Error: Debes especificar el paquete. Ejemplo: make upload pkg=QMB)
endif
	@VER=$$(cat $(pkg)/version.txt); \
	Filename="$(pkg)_v$${VER}_Source.txt"; \
	if [ ! -f "$$Filename" ]; then \
		echo "Archivo $$Filename no encontrado. Generándolo..."; \
		$(MAKE) $(pkg); \
	fi; \
	echo "Subiendo $$Filename a Google Drive ($(DrivePath))..."; \
	rclone copy "$$Filename" $(DriveRemote):$(DrivePath) --progress && \
	echo "Carga de $$Filename completada con éxito."

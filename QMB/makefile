# Variables de configuración (CamelCase obligatorio)
OutputSource = QMB_Source.txt
FileExtensions = -name "*.m" -o -name "*.wl"

# Objetos virtuales
.PHONY: all clean $(OutputSource)

# Objetivo principal
all: $(OutputSource)

$(OutputSource):
	@echo "Iniciando la extracción del código fuente de QMB..."
	@find . -type f \( $(FileExtensions) \) | while read -r File; do \
		echo "=== FILENAME: QMB/$${File#./} ==="; \
		cat "$$File"; \
		echo "-e "; \
		echo ""; \
	done > $(OutputSource)
	@echo "¡Archivo $(OutputSource) generado exitosamente en el directorio!"

clean:
	@echo "Eliminando $(OutputSource)..."
	rm -f $(OutputSource)

# Uso: 
# make bump type=major
# make bump type=minor
# make bump type=patch

bump:
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
		echo "Error: Define type como 'major', 'minor' o 'patch'"; \
		exit 1; \
	fi; \
	NewVersion="$$Major.$$Minor.$$Patch"; \
	echo $$NewVersion > version.txt; \
	sed -i "s/^\([ \t]*\)Version -> \"[^\"]*\"/\1Version -> \"$$NewVersion\"/" PacletInfo.wl; \
	echo "Actualización $(type) exitosa. Nueva versión: $$NewVersion"

# Exit with nonzero exit code if anything fails
set -e

# Remove the current documentation
rm -rf documentation/html

# Generate the Doxygen code documentation and log the output.
echo 'Generating Doxygen code documentation...'

# Redirect both stderr and stdout to the log file AND the console.
doxygen documentation/spipack.doxyfile.in 2>&1 | tee doxygen.log

# Upload the documentation if Doxygen successfully created the documentation.
if [ -d "documentation/html" ] && [ -f "documentation/html/index.html" ]; then
    echo 'Uploading documentation to the documentation branch...'
    # Add everything in the documentation directory
    git add documentation/html

    # Commit the added files with a title and description containing the Travis CI
    git commit -m "Deploy code docs to GitHub Pages Travis build: ${TRAVIS_BUILD_NUMBER} [ci skip]"

    # Push to the remote master branch.
    git push https://${GH_REPO_TOKEN}@github.com/SuperParameterIzation/spipack.git HEAD:master
else
    echo '' >&2
    echo 'Warning: No documentation (html) files have been found!' >&2
    echo 'Warning: Not going to push the documentation to GitHub!' >&2
    exit 1
fi

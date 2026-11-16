# Keep only the intended top-level HTML pages in a generated pkgdown site.
# Run after building the site: Rscript tools/filter-pkgdown-pages.R [site directory]
# Reference and article pages in subdirectories are not filtered.
args <- commandArgs(trailingOnly = TRUE)
site <- if (length(args)) args[[1]] else "docs"
if (!dir.exists(site)) stop("Site directory does not exist: ", site)

allowed <- c("index.html", "authors.html", "LICENSE.html", "404.html")
pages <- list.files(site, pattern = "\\.html$", full.names = TRUE, ignore.case = TRUE)
unwanted <- pages[!basename(pages) %in% allowed]
if (length(unwanted)) {
  message("Removing unwanted top-level pages: ", paste(basename(unwanted), collapse = ", "))
  if (any(unlink(unwanted) != 0L)) stop("Could not remove unwanted pages")
  companions <- file.path(site, sub("\\.html$", ".md", basename(unwanted), ignore.case = TRUE))
  if (any(unlink(companions[file.exists(companions)]) != 0L))
    stop("Could not remove unwanted companion pages")
}

# Refresh the search index so it cannot retain content from removed pages.
if (file.exists(file.path(site, "search.json")))
  pkgdown::build_search(override = list(destination = normalizePath(site)))

# Remove references to deleted pages from the existing sitemap.
sitemap <- file.path(site, "sitemap.xml")
if (length(unwanted) && file.exists(sitemap)) {
  doc <- xml2::read_xml(sitemap)
  entries <- xml2::xml_find_all(doc, ".//*[local-name()='url']")
  locations <- xml2::xml_text(xml2::xml_find_first(entries, "./*[local-name()='loc']"))
  remove <- vapply(locations, function(url) {
    any(endsWith(url, paste0("/", basename(unwanted))))
  }, logical(1))
  xml2::xml_remove(entries[remove])
  xml2::write_xml(doc, sitemap)
}

document.addEventListener('DOMContentLoaded', () => {
    const navLinks = document.querySelectorAll('a.md-nav__link , a.md-tabs__link');
    navLinks.forEach(link => {
        link.addEventListener('click', event => {
            // Force a full page reload
            window.location.href = event.currentTarget.href;
            event.preventDefault();
        });
    });
    // Handle browser back/forward navigation
    window.addEventListener('popstate', () => {
        // Reload the page when navigating using the Back or Forward buttons
        window.location.reload();
    });
});
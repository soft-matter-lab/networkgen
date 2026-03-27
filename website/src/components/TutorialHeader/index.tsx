import React from 'react';
import styles from './styles.module.css';

type Tag = {
  label: string;
};

type TutorialHeaderProps = {
  title: string;
  description: string;
  imageSrc?: string;
  imagePlaceholder?: string;
  tags?: Tag[];
  difficulty?: 'beginner' | 'intermediate' | 'advanced';
};

const difficultyColor = {
  beginner:     { bg: '#E1F5EE', color: '#085041' },
  intermediate: { bg: '#FAEEDA', color: '#633806' },
  advanced:     { bg: '#FAECE7', color: '#712B13' },
};

export default function TutorialHeader({
  title,
  description,
  imageSrc,
  imagePlaceholder = 'Network image',
  tags = [],
  difficulty,
}: TutorialHeaderProps) {
  const dc = difficulty ? difficultyColor[difficulty] : null;
  return (
    <div className={styles.header}>
      <div className={styles.imageCol}>
        {imageSrc
          ? <img src={imageSrc} alt={title} className={styles.image} />
          : <div className={styles.imagePlaceholder}>{imagePlaceholder}</div>
        }
      </div>
      <div className={styles.metaCol}>
        <p className={styles.metaLabel}>What you'll build</p>
        <p className={styles.metaDesc}>{description}</p>
        <div className={styles.tags}>
          {difficulty && dc && (
            <span className={styles.tag} style={{ background: dc.bg, color: dc.color }}>
              {difficulty}
            </span>
          )}
          {tags.map(t => (
            <span key={t.label} className={styles.tag}>{t.label}</span>
          ))}
        </div>
      </div>
    </div>
  );
}
